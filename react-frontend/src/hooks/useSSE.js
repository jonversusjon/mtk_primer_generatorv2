// useSSE.js
import { useEffect, useState, useRef } from "react";
import { SSE_BASE_URL } from "../config/config.js";

const useSSE = (jobId, sequenceIdx) => {
  const [sseEvent, setSseEvent] = useState(null);
  const eventSourceRef = useRef(null);
  const reconnectTimeoutRef = useRef(null);
  const connectAttemptRef = useRef(0);

  // Build the channel URL. The channel should match the naming convention used by your backend.
  const channel = `${SSE_BASE_URL}/stream?channel=job_${jobId}_${sequenceIdx}`;

  useEffect(() => {
    if (!jobId || sequenceIdx === undefined) {
      console.warn("SSE Hook: Missing jobId or sequenceIdx, not connecting");
      return;
    }

    // Clear any pending reconnect attempts
    if (reconnectTimeoutRef.current) {
      clearTimeout(reconnectTimeoutRef.current);
    }

    const connect = () => {
      console.log(
        `SSE Hook: Connecting to ${channel} (attempt ${connectAttemptRef.current})`
      );

      // Close any existing connection
      if (eventSourceRef.current) {
        eventSourceRef.current.close();
      }

      const eventSource = new EventSource(channel);
      eventSourceRef.current = eventSource;

      eventSource.onopen = () => {
        console.log(`SSE Hook: Connection opened to ${channel}`);
        connectAttemptRef.current = 0; // Reset attempt counter
      };

      eventSource.onmessage = (event) => {
        try {
          const data = JSON.parse(event.data);
          // Allow data to be either wrapped inside a "data" property or as a plain object
          const parsed = (typeof data === "object" && data.data) ? data.data : data;
          if (!parsed.timestamp) {
            parsed.timestamp = Date.now();
          }
          console.log(`SSE Hook: Received data from ${channel}:`, parsed);
          setSseEvent(parsed);
        } catch (error) {
          console.error(
            `SSE Hook: Error parsing SSE data from ${channel}:`,
            error,
            event.data
          );
        }
      };

      eventSource.onerror = (error) => {
        console.error(`SSE Hook: Error on ${channel}:`, error);
        eventSource.close();

        // Exponential backoff for reconnection; max delay 30 seconds
        const delay = Math.min(1000 * 2 ** connectAttemptRef.current, 30000);
        console.log(`SSE Hook: Reconnecting in ${delay}ms...`);
        connectAttemptRef.current++;
        reconnectTimeoutRef.current = setTimeout(() => {
          connect();
        }, delay);
      };
    };

    connect();

    // Cleanup on unmount
    return () => {
      console.log(`SSE Hook: Cleanup - closing connection to ${channel}`);
      if (eventSourceRef.current) {
        eventSourceRef.current.close();
      }
      if (reconnectTimeoutRef.current) {
        clearTimeout(reconnectTimeoutRef.current);
      }
    };
  }, [jobId, sequenceIdx, channel]);

  return sseEvent;
};

export default useSSE;
