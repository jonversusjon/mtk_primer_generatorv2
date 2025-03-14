import { API_BASE_URL } from "../config/config.js";

/**
 * Fetch available species for the dropdown
 * @returns {Promise<Array>} Array of available species
 */
export const fetchAvailableSpecies = async () => {
  try {
    const response = await fetchWithErrorHandling("/api/species");
    console.log("Available species:", response);
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
  console.group("Protocol Generation Process");
  console.log(
    "Starting protocol generation with data:",
    JSON.stringify(formData, null, 2)
  );

  try {
    console.log(`Sending request to ${API_BASE_URL}/generate_protocol`);
    console.time("Protocol generation request");

    const response = await fetch(`${API_BASE_URL}/generate_protocol`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify(formData),
    });

    console.timeEnd("Protocol generation request");
    console.log(`Response status: ${response.status} ${response.statusText}`);

    // Log response headers
    const headers = {};
    response.headers.forEach((value, key) => {
      headers[key] = value;
    });
    console.log("Response headers:", headers);

    console.log("Processing JSON response");
    const data = await response.json();
    console.log("Response data:", data);

    if (!response.ok) {
      console.error("Server returned error:", data);
      throw new Error(data.error || "Failed to generate protocol");
    }

    console.log("Protocol generation successful");
    console.groupEnd();
    return data;
  } catch (error) {
    console.error("Protocol generation failed:", error);
    console.error("Error stack:", error.stack);
    console.groupEnd();
    throw error;
  }
};

/**
 * Export protocol results as a TSV file
 * @param {Object} protocolData - Protocol data to export
 * @returns {Promise<Object>} Response with download URL
 */
export const exportProtocolAsTsv = async (protocolData) => {
  console.group("Protocol Export Process");
  console.log(
    "Exporting protocol data:",
    JSON.stringify(protocolData, null, 2)
  );

  try {
    console.time("Protocol export request");
    const data = await fetchWithErrorHandling("/export_protocol", {
      method: "POST",
      body: JSON.stringify(protocolData),
    });
    console.timeEnd("Protocol export request");

    console.log("Export successful, response:", data);
    console.groupEnd();
    return data;
  } catch (error) {
    console.error("Protocol export failed:", error);
    console.error("Error stack:", error.stack);
    console.groupEnd();
    throw error;
  }
};

// Utility function to add logging to the fetchWithErrorHandling function
// Note: If fetchWithErrorHandling is defined elsewhere, you may want to modify it directly instead
const fetchWithErrorHandling = async (url, options = {}) => {
  console.group(`API Request: ${options.method || "GET"} ${url}`);
  console.log("Request options:", options);

  try {
    console.time("Request execution");
    const response = await fetch(url, options);
    console.timeEnd("Request execution");

    console.log(`Response status: ${response.status} ${response.statusText}`);

    // Log headers
    const headers = {};
    response.headers.forEach((value, key) => {
      headers[key] = value;
    });
    console.log("Response headers:", headers);

    // Clone the response for logging (since response body can only be read once)
    const responseClone = response.clone();

    // Process response based on content type
    const contentType = response.headers.get("content-type");
    console.log(`Response content type: ${contentType}`);

    if (contentType && contentType.includes("application/json")) {
      const responseData = await response.json();
      console.log("Response data:", responseData);

      if (!response.ok) {
        throw new Error(
          responseData.error || `Request failed with status ${response.status}`
        );
      }

      console.groupEnd();
      return responseData;
    } else {
      const text = await responseClone.text();
      console.log(`Non-JSON response (${text.length} chars)`);

      if (!response.ok) {
        throw new Error(`Request failed with status ${response.status}`);
      }

      console.groupEnd();
      return text;
    }
  } catch (error) {
    console.error("Request failed:", error);
    console.error("Error stack:", error.stack);
    console.groupEnd();
    throw error;
  }
};
