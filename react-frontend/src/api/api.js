import { API_BASE_URL } from "../config/config.js";

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
