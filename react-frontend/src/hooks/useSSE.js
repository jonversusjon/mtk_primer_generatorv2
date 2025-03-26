import { useEffect, useState } from "react";
import { SSE_BASE_URL } from "../config/config.js";

const useSSE = (jobId, sequenceIdx) => {
  const [sseEvent, setSseEvent] = useState(null);
  // Make sure the channel name matches exactly what's used in celery_tasks.py
  const channel = `${SSE_BASE_URL}/stream?channel=job_${jobId}_${sequenceIdx}`;

  useEffect(() => {
    if (!jobId || sequenceIdx === undefined) {
      console.warn("SSE Hook: Missing jobId or sequenceIdx, not connecting");
      return;
    }

    console.log(`SSE Hook: Connecting to channel ${channel}`);
    const eventSource = new EventSource(channel);

    eventSource.onmessage = (event) => {
      try {
        // Check the response structure
        const eventData = JSON.parse(event.data);
        console.log("Raw SSE data received:", eventData);
        
        // Handle different possible structures
        // Sometimes Flask-SSE wraps the data in a data property
        const parsed = eventData.data || eventData;
        
        console.log("Processed SSE in hook on channel", channel, ":", parsed);
        setSseEvent(parsed);
      } catch (error) {
        console.error("Error parsing SSE data:", error, event.data);
      }
    };

    eventSource.onerror = (error) => {
      console.error("SSE error on channel", channel, ":", error);
      // Implement reconnection with exponential backoff
      // This helps ensure we don't miss events due to connection issues
    };

    // Cleanup: only close when the component using this hook unmounts
    return () => {
      console.log("Closing SSE connection for job", jobId);
      eventSource.close();
    };
  }, [jobId, sequenceIdx, channel]);

  return sseEvent;
};

export default useSSE;