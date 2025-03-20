import { API_BASE_URL } from "../config/config.js";

// Helper function to delay for a specified number of milliseconds.
const delay = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

/**
 * Constructs an SSE URL based on the API_BASE_URL without adding an extra /api.
 * @param {string} path - The SSE endpoint path (e.g., '/status/jobId').
 * @returns {string} - The full URL for SSE connection.
 */
const constructSSEUrl = (path) => {
  // Remove trailing slash if present.
  const base = API_BASE_URL.replace(/\/$/, "");
  // Directly append the path to the base. This avoids any extra /api prefix.
  return `${base}${path}`;
};

/**
 * Utility function to make API requests with error handling and detailed logging.
 * @param {string} url - The endpoint URL (relative or absolute).
 * @param {Object} options - Fetch options (method, headers, body, etc.).
 * @returns {Promise<Object|string>} The parsed JSON response or text.
 */
const fetchWithErrorHandling = async (url, options = {}) => {
  const fullUrl = url.startsWith("http") ? url : `${API_BASE_URL}${url}`;
  console.group(`API Request: ${options.method || "GET"} ${fullUrl}`);
  console.log("Request options:", options);

  try {
    console.time("Request execution");
    const response = await fetch(fullUrl, options);
    console.timeEnd("Request execution");
    console.log(`Response status: ${response.status} ${response.statusText}`);

    // Log response headers.
    const headers = {};
    response.headers.forEach((value, key) => (headers[key] = value));
    console.log("Response headers:", headers);

    // Process response based on content type.
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
      const text = await response.text();
      console.log(`Non-JSON response (${text.length} chars)`);
      if (!response.ok) {
        throw new Error(`Request failed with status ${response.status}`);
      }
      console.groupEnd();
      return text;
    }
  } catch (error) {
    console.error("Request failed:", error);
    console.groupEnd();
    throw error;
  }
};

/**
 * Fetch available species for the dropdown.
 * @returns {Promise<Array>} Array of available species.
 */
export const fetchAvailableSpecies = async () => {
  try {
    const response = await fetchWithErrorHandling("/species");
    console.log("Available species:", response);
    return response;
  } catch (error) {
    console.error("Error fetching species:", error);
    throw error;
  }
};

/**
 * Validate a DNA sequence.
 * @param {string} sequence - DNA sequence to validate.
 * @returns {Promise<Object>} Validation result.
 */
export const validateSequence = async (sequence) => {
  try {
    const data = await fetchWithErrorHandling("/validation/validate/sequence", {
      method: "POST",
      body: JSON.stringify({ sequence }),
    });
    return data;
  } catch (error) {
    console.error("Error validating sequence:", error);
    throw error;
  }
};

/**
 * Establish an SSE connection to monitor protocol generation progress.
 * Supports auto-reconnection on errors with exponential backoff.
 * @param {string} jobId - The job ID to monitor.
 * @param {Function} onStatusUpdate - Callback for status updates.
 * @param {Object} options - Optional settings (e.g., maxRetries, baseDelay).
 * @returns {EventSource} The event source object.
 */
export const monitorProtocolProgress = (jobId, onStatusUpdate, options = {}) => {
  const maxRetries = options.maxRetries || 5;
  const baseDelay = options.baseDelay || 2000; // Base delay in milliseconds
  let retries = 0;
  let eventSource = null;

  const connect = () => {
    console.log(
      `Setting up SSE connection for job ${jobId} (Attempt ${retries + 1})`
    );
    // Build SSE URL without adding an extra /api prefix.
    const sseUrl = constructSSEUrl(`/status/${jobId}`);
    eventSource = new EventSource(sseUrl);

    eventSource.onmessage = (event) => {
      try {
        const statusData = JSON.parse(event.data);
        console.log(`Progress update for job ${jobId}:`, statusData);
        onStatusUpdate(statusData);

        // Close connection if job is complete or an error occurred.
        if (statusData.percentage === 100 || statusData.percentage === -1) {
          console.log(`Job ${jobId} complete or failed, closing SSE connection`);
          eventSource.close();
        }
      } catch (error) {
        console.error("Error processing SSE message:", error);
      }
    };

    eventSource.onerror = (error) => {
      console.error("SSE connection error:", error);
      eventSource.close();
      if (retries < maxRetries) {
        retries++;
        // Exponential backoff: delay = baseDelay * 2^(retries - 1)
        const reconnectDelay = baseDelay * Math.pow(2, retries - 1);
        console.log(`Attempting to reconnect in ${reconnectDelay}ms...`);
        setTimeout(() => {
          connect();
        }, reconnectDelay);
      } else {
        console.error("Max retries reached for SSE connection");
        onStatusUpdate({
          message: "Connection to server lost. The job is still processing.",
          percentage: -1,
          step: "error",
        });
      }
    };

    return eventSource;
  };

  return connect();
};

/**
 * Generate a protocol based on form data.
 * Optionally sets up an SSE connection for status updates.
 * Falls back to polling if no onStatusUpdate callback is provided.
 * @param {Object} formData - Form data containing sequence info.
 * @param {Function} onStatusUpdate - Optional callback for status updates.
 * @returns {Promise<Object>} Protocol results or initial data with eventSource.
 */
export const generateProtocol = async (formData, onStatusUpdate) => {
  console.group("Protocol Generation Process");
  console.log(
    "Starting protocol generation with data:",
    JSON.stringify(formData, null, 2)
  );

  // Generate a job ID if not provided.
  const jobId = formData.jobId || Date.now().toString();
  const dataWithJobId = { ...formData, jobId };
  let eventSource = null;

  try {
    // Set up SSE connection if a status callback is provided.
    if (onStatusUpdate) {
      eventSource = monitorProtocolProgress(jobId, onStatusUpdate);
    }

    console.log(`Sending request to ${API_BASE_URL}/generate_protocol`);
    console.time("Protocol generation request");

    const response = await fetch(`${API_BASE_URL}/generate_protocol`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(dataWithJobId),
    });

    console.timeEnd("Protocol generation request");
    console.log(`Response status: ${response.status} ${response.statusText}`);

    // Log response headers.
    const headers = {};
    response.headers.forEach((value, key) => (headers[key] = value));
    console.log("Response headers:", headers);

    console.log("Processing JSON response");
    const initialData = await response.json();
    console.log("Initial response data:", initialData);

    if (!response.ok) {
      console.error("Server returned error:", initialData);
      throw new Error(initialData.error || "Failed to generate protocol");
    }

    // Store submitted data for immediate UI rendering.
    sessionStorage.setItem(
      "results",
      JSON.stringify(initialData.submittedData)
    );

    // Fall back to polling if no onStatusUpdate callback is provided.
    if (initialData.jobId && !onStatusUpdate) {
      console.log(`Polling for results of job ${initialData.jobId}`);
      const finalData = await pollForResults(initialData.jobId);
      console.log("Final data from polling:", finalData);
      console.groupEnd();
      return finalData;
    }

    console.log("Protocol generation initiated successfully");
    console.groupEnd();
    return { initialData, eventSource };
  } catch (error) {
    console.error("Protocol generation failed:", error);
    console.groupEnd();
    throw error;
  }
};

/**
 * Poll for final results of a protocol generation job.
 * @param {string} jobId - The job ID to poll for.
 * @param {number} maxAttempts - Maximum polling attempts.
 * @param {number} interval - Interval between polls (ms).
 * @returns {Promise<Object>} Final protocol results.
 */
export const pollForResults = async (
  jobId,
  maxAttempts = 30,
  interval = 2000
) => {
  console.log(`Starting to poll for results of job ${jobId}`);

  for (let attempt = 1; attempt <= maxAttempts; attempt++) {
    console.log(`Poll attempt ${attempt}/${maxAttempts}`);

    try {
      const response = await fetch(`${API_BASE_URL}/results/${jobId}`);
      const data = await response.json();
      console.log("Poll result:", data);

      if (response.ok && data.complete) {
        console.log(`Job ${jobId} complete, returning results`);
        return data.data;
      }
      await delay(interval);
    } catch (error) {
      console.error(`Error polling for results on attempt ${attempt}:`, error);
    }
  }

  throw new Error(
    `Timed out waiting for results after ${maxAttempts} attempts`
  );
};

/**
 * Get the current status of a job.
 * @param {string} jobId - The job ID to check.
 * @returns {Promise<Object>} Current job status.
 */
export const getJobStatus = async (jobId) => {
  try {
    const response = await fetch(`${API_BASE_URL}/results/${jobId}`);
    const data = await response.json();
    return data;
  } catch (error) {
    console.error("Error getting job status:", error);
    throw error;
  }
};

/**
 * Export protocol results as a TSV file.
 * @param {Object} protocolData - Protocol data to export.
 * @returns {Promise<Object>} Response with the download URL.
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
    console.groupEnd();
    throw error;
  }
};
