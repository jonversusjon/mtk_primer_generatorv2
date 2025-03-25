import { useEffect, useState } from "react";
import { SSE_BASE_URL } from "../config/config.js";

const useSSE = (jobId, sequenceIdx) => {
  const [sseEvent, setSseEvent] = useState(null);

  useEffect(() => {
    const eventSource = new EventSource(
      `${SSE_BASE_URL}/stream?channel=job_${jobId}_${sequenceIdx}`
    );

    eventSource.onmessage = (event) => {
      console.log("Received SSE event:", event);
      try {
        const parsedData = JSON.parse(event.data);
        setSseEvent(parsedData); // update with the latest event
      } catch (err) {
        console.error("Error parsing SSE data:", err);
      }
    };

    eventSource.onerror = (error) => {
      console.error("SSE error:", error);
      eventSource.close();
    };

    return () => {
      console.log("Closing SSE connection");
      eventSource.close();
    };
  }, [jobId, sequenceIdx]);

  return sseEvent;
};

export default useSSE;
