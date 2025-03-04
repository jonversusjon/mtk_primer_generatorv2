import { API_BASE_URL } from "../config/config.js";

/**
 * API functions for communicating with the Flask backend
 */

/**
 * Helper function to handle fetch requests
 * @param {string} endpoint - API endpoint to call
 * @param {Object} options - Fetch options (method, headers, body)
 * @returns {Promise} Promise that resolves to the JSON response
 */
export const fetchWithErrorHandling = async (url, options = {}) => {
  try {
    const response = await fetch(url, options);
    if (!response.ok) {
      throw new Error(`API request failed with status ${response.status}`);
    }
    return await response.json();
  } catch (error) {
    console.error("API request failed:", error);
    throw error;
  }
};

/**
 * Fetch available species for the dropdown
 * @returns {Promise<Array>} Array of available species
 */
export const fetchAvailableSpecies = async () => {
  try {
    const response = await fetchWithErrorHandling("/api/species");
    return response;
  } catch (error) {
    console.error("Error fetching species:", error);
    throw error;
  }
};

/**
 * Validate a DNA sequence
 * @param {string} sequence - DNA sequence to validate
 * @returns {Promise<Object>} Validation result
 */
export const validateSequence = async (sequence) => {
  const data = await fetchWithErrorHandling("/validation/validate/sequence", {
    method: "POST",
    body: JSON.stringify({ sequence }),
  });

  return data;
};

/**
 * Generate a protocol based on form data
 * @param {Object} formData - Form data containing sequence info
 * @returns {Promise<Object>} Protocol results
 */
export const generateProtocol = async (formData) => {
  // Create FormData object for file uploads
  const form = new FormData();

  // Convert sequences array to the format expected by the server
  formData.sequences.forEach((seq, index) => {
    form.append(`sequences[${index}][sequence]`, seq.sequence);
    form.append(`sequences[${index}][primerName]`, seq.primerName);
    form.append(`sequences[${index}][mtkPart]`, seq.mtkPart);
  });

  // Add other form fields
  form.append("numSequences", formData.numSequences);
  form.append("templateSequence", formData.templateSequence);
  form.append("species", formData.species);
  form.append("kozak", formData.kozak);
  form.append("max_mut_per_site", formData.max_mut_per_site);

  if (formData.verbose_mode) {
    form.append("verbose_mode", "on");
  }

  try {
    const response = await fetch(`${API_BASE_URL}/generate_protocol`, {
      method: "POST",
      body: form,
    });

    // Check if the response is JSON
    const contentType = response.headers.get("content-type");
    if (contentType && contentType.includes("application/json")) {
      const data = await response.json();

      if (!response.ok) {
        throw new Error(data.error || "Failed to generate protocol");
      }

      return data;
    } else {
      // Handle HTML response (error page)
      const text = await response.text();

      // Extract error message from HTML if possible
      const errorMatch = text.match(
        /<div class="alert alert-danger">(.*?)<\/div>/
      );
      if (errorMatch && errorMatch[1]) {
        throw new Error(errorMatch[1].trim());
      } else {
        throw new Error("Server returned an unexpected response");
      }
    }
  } catch (error) {
    console.error("Protocol generation failed:", error);
    throw error;
  }
};

/**
 * Export protocol results as a TSV file
 * @param {Object} protocolData - Protocol data to export
 * @returns {Promise<Object>} Response with download URL
 */
export const exportProtocolAsTsv = async (protocolData) => {
  const data = await fetchWithErrorHandling("/export_protocol", {
    method: "POST",
    body: JSON.stringify(protocolData),
  });

  return data;
};
