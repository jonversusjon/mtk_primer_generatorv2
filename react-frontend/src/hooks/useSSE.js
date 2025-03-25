import { useEffect } from "react";
import { SSE_BASE_URL } from "../config/config.js";

const useSSE = (jobId, sequenceIdx) => {
  useEffect(() => {
    const eventSource = new EventSource(
      `${SSE_BASE_URL}/stream?channel=job_${jobId}_${sequenceIdx}`
    );

    eventSource.onmessage = (event) => {
      console.log("Received SSE event:", event.data);
      // Process the event data as needed
    };

    eventSource.onerror = (error) => {
      console.error("SSE error:", error);
      eventSource.close();
    };

    return () => {
      eventSource.close();
    };
  }, [jobId, sequenceIdx]);
};

export default useSSE;
