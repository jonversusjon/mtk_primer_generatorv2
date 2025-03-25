import { useEffect, useState } from "react";
import { SSE_BASE_URL } from "../config/config.js";

const useSSE = (jobId, sequenceIdx) => {
  const [sseData, setSseData] = useState(null);

  useEffect(() => {
    const eventSource = new EventSource(
      `${SSE_BASE_URL}/stream?channel=job_${jobId}_${sequenceIdx}`
    );

    eventSource.onmessage = (event) => {
      console.log("Received SSE event:", event.data);
      try {
        const parsedData = JSON.parse(event.data);
        setSseData(parsedData);
      } catch (err) {
        console.error("Error parsing SSE data:", err);
      }
    };

    eventSource.onerror = (error) => {
      console.error("SSE error:", error);
      eventSource.close();
    };

    return () => {
      eventSource.close();
    };
  }, [jobId, sequenceIdx]);

  return sseData;
};

export default useSSE;
