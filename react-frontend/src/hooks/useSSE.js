// hooks/useSSE.js
// ... (imports and setup remain the same) ...
import { useEffect, useRef } from "react";
import { SSE_BASE_URL } from "../config/config.js";

const useSSE = (jobId, sequenceIdx, onMessage) => {
  const eventSourceRef = useRef(null);
  const reconnectTimeoutRef = useRef(null);
  const connectAttemptRef = useRef(0);
  const instanceIdRef = useRef(Math.random().toString(36).substring(7));
  const onMessageRef = useRef(onMessage);

  useEffect(() => {
    onMessageRef.current = onMessage;
  }, [onMessage]);

  const channel = `${SSE_BASE_URL}/stream?channel=job_${jobId}_${sequenceIdx}`;
  const logPrefix = `SSE Hook [${instanceIdRef.current} - ${sequenceIdx}]`;

  useEffect(() => {
    console.log(`${logPrefix}: useEffect RUNNING. JobId: ${jobId}, SeqIdx: ${sequenceIdx}`);
    // ... (checks for jobId, sequenceIdx, onMessageRef remain the same) ...
    if (!jobId || sequenceIdx === undefined) return;
    if (typeof onMessageRef.current !== 'function') return;


    let currentEventSource = null;

    const connect = () => {
      // ... (connect logic remains the same) ...
        if (eventSourceRef.current && eventSourceRef.current.readyState !== EventSource.CLOSED) {
            console.log(`${logPrefix}: connect() called, but EventSource exists and is not closed (State: ${eventSourceRef.current.readyState}). Aborting new connection.`);
            return;
        }
        console.log(`${logPrefix}: Connecting to ${channel} (attempt ${connectAttemptRef.current})`);
        const eventSource = new EventSource(channel);
        currentEventSource = eventSource;
        eventSourceRef.current = eventSource;


      eventSource.onopen = () => {
        // ... (onopen logic remains the same) ...
        console.log(
          `${logPrefix}: **Connection opened**. State: ${eventSource.readyState}`
        );
        connectAttemptRef.current = 0;
        if (reconnectTimeoutRef.current) {
          clearTimeout(reconnectTimeoutRef.current);
          reconnectTimeoutRef.current = null;
          console.log(
            `${logPrefix}: Cleared pending reconnect on connection open.`
          );
        }
      };

      eventSource.onmessage = (event) => {
         console.log(`${logPrefix}: Message received (raw). State: ${eventSource.readyState}`);
        try {
          const data = JSON.parse(event.data);
          const parsed =
            typeof data === "object" &&
            data !== null &&
            data.hasOwnProperty("data")
              ? data.data
              : data;

          if (
            typeof parsed === "object" &&
            parsed !== null &&
            !parsed.timestamp
          ) {
            parsed.timestamp = Date.now();
          }

           // *** ADDED DETAILED LOGGING AROUND CALLBACK ***
           console.log(`${logPrefix}: Parsed data in onmessage:`, parsed);
           console.log(`${logPrefix}: Type of onMessageRef.current: ${typeof onMessageRef.current}`);
           if (onMessageRef.current) {
               console.log(`${logPrefix}: Attempting to call onMessageRef.current...`);
               onMessageRef.current(parsed); // <-- The call
               console.log(`${logPrefix}: Successfully called onMessageRef.current.`); // Log AFTER the call
           } else {
               console.error(`${logPrefix}: onMessageRef.current is null or undefined! Cannot call callback.`);
           }
           // *** END ADDED LOGGING ***

        } catch (error) {
          console.error(`${logPrefix}: Error parsing JSON:`, error, event.data);
        }
      };

      eventSource.onerror = (error) => {
         // ... (onerror logic remains the same) ...
          console.error(
           `${logPrefix}: **onerror triggered**. State: ${eventSource.readyState}. Error object:`, error
          );
          if (error instanceof Event) {
             console.error(
               `${logPrefix}: onerror Event details: type=${error.type}, targetState=${error.target?.readyState}`
             );
           }
           if (eventSource.readyState !== EventSource.CLOSED) {
             console.log(`${logPrefix}: Closing EventSource due to error.`);
             eventSource.close();
           }
           if (eventSourceRef.current === eventSource) {
             eventSourceRef.current = null;
           }
           if (reconnectTimeoutRef.current) {
             clearTimeout(reconnectTimeoutRef.current);
             reconnectTimeoutRef.current = null;
           }
           const delay = Math.min(1000 * 2 ** connectAttemptRef.current, 30000);
           console.log(`${logPrefix}: Scheduling reconnect via onerror in ${delay}ms...`);
           connectAttemptRef.current++;
           reconnectTimeoutRef.current = setTimeout(connect, delay);
      };
    };

    // ... (Initial connect call and cleanup logic remain the same) ...
    console.log(`${logPrefix}: Calling connect(). Current source state: ${eventSourceRef.current?.readyState ?? 'N/A'}`);
    connect();

    return () => {
       console.log(
        `${logPrefix}: useEffect CLEANUP RUNNING. Closing connection if open. Source state: ${currentEventSource?.readyState ?? 'N/A'}`
       );
       // ... rest of cleanup
        if (
            currentEventSource &&
            currentEventSource.readyState !== EventSource.CLOSED
        ) {
            console.log(`${logPrefix}: Cleanup closing EventSource instance.`);
            currentEventSource.close();
        }
        if (reconnectTimeoutRef.current) {
            console.log(`${logPrefix}: Cleanup clearing reconnect timeout.`);
            clearTimeout(reconnectTimeoutRef.current);
            reconnectTimeoutRef.current = null;
        }
        if (eventSourceRef.current === currentEventSource) {
            console.log(`${logPrefix}: Cleanup nullifying eventSourceRef.`);
            eventSourceRef.current = null;
        }
    };
  }, [jobId, sequenceIdx, channel, logPrefix]); // Keep onMessage removed

  return null;
};

export default useSSE;