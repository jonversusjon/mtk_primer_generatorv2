import { useEffect, useState, useRef } from "react";
import { SSE_BASE_URL } from "../config/config.js";

const useSSE = (jobId, sequenceIdx) => {
  const [sseEvent, setSseEvent] = useState(null);
  const eventSourceRef = useRef(null);
  const reconnectTimeoutRef = useRef(null);
  const connectAttemptRef = useRef(0);

  // Make sure the channel name matches exactly what's used in celery_tasks.py
  const channel = `${SSE_BASE_URL}/stream?channel=job_${jobId}_${sequenceIdx}`;

  useEffect(() => {
    if (!jobId || sequenceIdx === undefined) {
      console.warn("SSE Hook: Missing jobId or sequenceIdx, not connecting");
      return;
    }

    // Clear any existing timeouts to prevent memory leaks
    if (reconnectTimeoutRef.current) {
      clearTimeout(reconnectTimeoutRef.current);
    }

    // Function to connect to SSE with exponential backoff
    const connect = () => {
      console.log(
        `SSE Hook: Connecting to channel ${channel} (attempt ${connectAttemptRef.current})`
      );

      // Close existing connection if any
      if (eventSourceRef.current) {
        eventSourceRef.current.close();
      }

      const eventSource = new EventSource(channel);
      eventSourceRef.current = eventSource;

      eventSource.onopen = () => {
        console.log(`SSE Hook: Connection opened to ${channel}`);
        connectAttemptRef.current = 0; // Reset attempt counter on successful connection
      };

      eventSource.onmessage = (event) => {
        try {
          // Check the response structure
          const eventData = JSON.parse(event.data);
          console.log(`SSE Hook: Raw data received on ${channel}:`, eventData);

          // Sometimes Flask-SSE wraps the data in a data property
          let parsed;
          if (typeof eventData === "object" && eventData !== null) {
            parsed = eventData.data || eventData;

            // Add timestamp if not present
            if (!parsed.timestamp) {
              parsed.timestamp = Date.now();
            }
          } else {
            parsed = { data: eventData, timestamp: Date.now() };
          }

          console.log(`SSE Hook: Processed data on ${channel}:`, parsed);
          setSseEvent(parsed);
        } catch (error) {
          console.error(
            `SSE Hook: Error parsing data on ${channel}:`,
            error,
            event.data
          );
        }
      };

      eventSource.onerror = (error) => {
        console.error(`SSE Hook: Error on ${channel}:`, error);

        // Close the current connection
        eventSource.close();

        // Implement exponential backoff for reconnection
        const delay = Math.min(1000 * 2 ** connectAttemptRef.current, 30000); // Max 30 second delay
        console.log(`SSE Hook: Reconnecting in ${delay}ms...`);

        connectAttemptRef.current++; // Increment attempt counter

        reconnectTimeoutRef.current = setTimeout(() => {
          connect(); // Try to reconnect
        }, delay);
      };
    };

    // Initial connection
    connect();

    // Cleanup function
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
