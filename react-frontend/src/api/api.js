import { API_BASE_URL } from "../config/config.js";

/**
 * Fetch wrapper with standardized logging + error handling.
 */
export const fetchWithErrorHandling = async (url, options = {}) => {
  const fullUrl = url.startsWith("http") ? url : `${API_BASE_URL}${url}`;
  console.group(`API Request: ${options.method || "GET"} ${fullUrl}`);
  console.log("Request options:", options);

  try {
    console.time("Request execution");
    const response = await fetch(fullUrl, options);
    console.timeEnd("Request execution");
    console.log(`Response status: ${response.status} ${response.statusText}`);

    const contentType = response.headers.get("content-type");
    const data = contentType?.includes("application/json")
      ? await response.json()
      : await response.text();

    if (!response.ok) {
      throw new Error(data.error || `Request failed (${response.status})`);
    }

    console.log("Response data:", data);
    console.groupEnd();
    return data;
  } catch (error) {
    console.error("Request failed:", error);
    console.groupEnd();
    throw error;
  }
};

/**
 * Submit protocol generation; returns jobId for SSE subscription.
 */
export const submitProtocol = async (formData) => {
  const jobId = formData.jobId || Date.now().toString();
  const payload = { ...formData, jobId };

  console.group("Submit Protocol");
  console.time("Protocol request");

  const response = await fetch(`${API_BASE_URL}/generate_protocol`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });

  console.timeEnd("Protocol request");
  const initialData = await response.json();

  if (!response.ok) {
    console.error("Protocol generation error:", initialData);
    throw new Error(initialData.error);
  }

  console.log("Protocol generation started, jobId:", jobId);
  console.groupEnd();
  return { jobId };
};

/** Remaining helpers — unchanged **/
export const fetchAvailableSpecies = async () => fetchWithErrorHandling("/species");
export const validateSequence = async (seq) =>
  fetchWithErrorHandling("/validation/validate/sequence", {
    method: "POST",
    body: JSON.stringify({ sequence: seq }),
  });
