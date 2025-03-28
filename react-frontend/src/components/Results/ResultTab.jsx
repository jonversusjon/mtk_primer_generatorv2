// Results/ResultTab.jsx
import React, { useState, useEffect, useCallback, useRef } from "react";
import useSSE from "../../hooks/useSSE";
import ProtocolTracker from "../../components/Results/ProtocolTracker";

const ResultTab = ({ jobId, sequenceIdx }) => {
  // --- State and Refs ---
  const [rawSseEvents, setRawSseEvents] = useState(() => {
    const savedEvents = sessionStorage.getItem(
      `sseEvents_${jobId}_${sequenceIdx}`
    );
    try {
      console.log(
        `[ResultTab:${sequenceIdx}] Initial load from sessionStorage: ${
          savedEvents ? savedEvents.length + " chars" : "null"
        }`
      );
      return savedEvents ? JSON.parse(savedEvents) : [];
    } catch (e) {
      console.error(
        `[ResultTab:${sequenceIdx}] Failed to parse saved SSE events:`,
        e
      );
      return [];
    }
  });

  // Protocol Tracker State
  const [protocolSteps, setProtocolSteps] = useState([]);
  const [protocolMessages, setProtocolMessages] = useState([]);
  const [protocolCallouts, setProtocolCallouts] = useState([]);
  const [protocolSseData, setProtocolSseData] = useState({});

  const processedEventIds = useRef(
    new Set(
      rawSseEvents.map(
        (event) =>
          `seq${sequenceIdx}-${JSON.stringify({
            ...event,
            _idTimestamp: event.clientTimestamp,
          })}`
      )
    )
  );
  const streamClosed = useRef(false);
  const hasReceivedData = useRef(rawSseEvents.length > 0);
  const isMounted = useRef(false);

  // --- Process SSE data into Protocol Tracker format ---
  useEffect(() => {
    console.log(
      `[ResultTab:${sequenceIdx}] Processing ${rawSseEvents.length} events for protocol data`
    );
    if (rawSseEvents.length > 0) {
      console.log("Sample event structure:", rawSseEvents[0]);
    }

    // Define expected protocol steps
    const expectedSteps = [
      "Preprocessing",
      "Restriction Sites",
      "Mutation Analysis",
      "Primer Design",
      "PCR Reaction Grouping",
    ];

    // Initialize steps with default values
    const steps = expectedSteps.map((name) => ({
      name,
      status: "active",
      stepProgress: 0,
      message: "",
      notificationCount: 0,
    }));

    const messages = [];
    const callouts = [];
    const sseData = {};

    // Initialize empty data structure for each step
    expectedSteps.forEach((step) => {
      sseData[step] = { timestamp: Date.now() };
    });

    // Process each event to update steps and sseData
    rawSseEvents.forEach((event) => {
      // Try to determine which step this event belongs to
      let stepName = null;

      // Check for explicit step property
      if (event.step) {
        stepName = event.step;
      }
      // Or infer from type and content
      else if (event.type) {
        if (event.type.includes("restriction") || event.sites) {
          stepName = "Restriction Sites";
        } else if (event.type.includes("mutation") || event.mutations) {
          stepName = "Mutation Analysis";
        } else if (
          event.type.includes("primer") ||
          event.primers ||
          event.edgePrimers ||
          event.mutPrimers
        ) {
          stepName = "Primer Design";
        } else if (
          event.type.includes("pcr") ||
          event.reactions ||
          event.pcrReactions
        ) {
          stepName = "PCR Reaction Grouping";
        }
      }

      // If we identified a step, update its status and data
      if (stepName) {
        // Find the step in our array
        const stepIndex = steps.findIndex((s) => s.name === stepName);
        if (stepIndex >= 0) {
          // Update step status based on event properties
          if (event.status) {
            steps[stepIndex].status = event.status;
          } else if (event.stepProgress > 0) {
            steps[stepIndex].status = "active";
          }

          // Update progress if provided
          if (event.stepProgress !== undefined) {
            steps[stepIndex].stepProgress = event.stepProgress;
          }

          // Update notification count if provided
          if (event.notificationCount !== undefined) {
            steps[stepIndex].notificationCount = event.notificationCount;
          }

          // Update message if provided
          if (event.message) {
            steps[stepIndex].message = event.message;
            messages.push({
              step: stepName,
              message: event.message,
              timestamp: event.timestamp || event.clientTimestamp,
            });
          }

          // Update callout if provided
          if (event.callout) {
            callouts[stepIndex] = event.callout;
            callouts.push({
              step: stepName,
              callout: event.callout,
            });
            steps[stepIndex].callout = event.callout;
          }
        }

        // Update the sseData for this step
        const stepData = sseData[stepName] || {};

        // Merge all properties from the event
        sseData[stepName] = {
          ...stepData,
          ...event, // Copy all event properties directly
          timestamp: event.timestamp || event.clientTimestamp || Date.now(),
        };

        // Handle specific data structures we know are used by ProtocolTracker
        if (event.sites) {
          sseData[stepName].sites = event.sites;
        }
        if (event.mutations) {
          sseData[stepName].mutations = event.mutations;
        }
        if (event.edgePrimers) {
          sseData[stepName].edgePrimers = event.edgePrimers;
        }
        if (event.mutPrimers) {
          sseData[stepName].mutPrimers = event.mutPrimers;
        }
        if (event.reactions || event.pcrReactions) {
          sseData[stepName].pcrReactions =
            event.reactions || event.pcrReactions;
        }
        if (event.callout) {
          sseData[stepName].callout = event.callout;
        }
        if (event.notifcation_count) {
          sseData[stepName].notificationCount = event.notification_count;
        }
      }
    });

    // Update protocol data state
    setProtocolSteps(steps);
    console.log(
      "STEPS DATA FOR PROTOCOL TRACKER:",
      JSON.stringify(steps, null, 2)
    );

    setProtocolMessages(messages);
    setProtocolCallouts(callouts);
    setProtocolSseData(sseData);

    console.log(`[ResultTab:${sequenceIdx}] Protocol data processed:`, {
      steps: steps.length,
      messages: messages.length,
      sseDataKeys: Object.keys(sseData),
    });
  }, [rawSseEvents, sequenceIdx]);

  // --- Callback for useSSE ---
  const processSseEvent = useCallback(
    (eventData) => {
      console.log(
        `[ResultTab:${sequenceIdx}] processSseEvent ENTERED. Data received:`,
        eventData
      );

      if (!eventData || typeof eventData !== "object") {
        console.warn(
          `[ResultTab:${sequenceIdx}] processSseEvent: EXITING - Invalid eventData type or null.`
        );
        return;
      }
      if (!hasReceivedData.current) {
        hasReceivedData.current = true;
        console.log(
          `[ResultTab:${sequenceIdx}] processSseEvent: Set hasReceivedData = true.`
        );
      }
      if (
        eventData.sequenceIdx !== undefined &&
        eventData.sequenceIdx !== sequenceIdx
      ) {
        console.warn(
          `[ResultTab:${sequenceIdx}] processSseEvent: Skipping event due to sequenceIdx mismatch. Event idx: ${eventData.sequenceIdx}`
        );
        return;
      }
      const eventIdSource = { ...eventData };
      const consistentTimestamp = eventData.timestamp || Date.now();
      eventIdSource._idTimestamp = consistentTimestamp;
      const eventId = `seq${sequenceIdx}-${JSON.stringify(eventIdSource)}`;
      if (processedEventIds.current.has(eventId)) {
        return;
      }
      processedEventIds.current.add(eventId);
      const eventToStore = {
        ...eventData,
        clientTimestamp: consistentTimestamp,
      };
      setRawSseEvents((prevEvents) => {
        if (!isMounted.current) {
          console.warn(
            `[ResultTab:${sequenceIdx}] processSseEvent: Attempted state update after unmount. Skipping.`
          );
          return prevEvents;
        }
        console.log(
          `[ResultTab:${sequenceIdx}] processSseEvent: Updating state with event timestamp: ${eventToStore.clientTimestamp}`
        );
        const updatedEvents = [...prevEvents, eventToStore];
        updatedEvents.sort(
          (a, b) => Number(a.clientTimestamp) - Number(b.clientTimestamp)
        );
        try {
          sessionStorage.setItem(
            `sseEvents_${jobId}_${sequenceIdx}`,
            JSON.stringify(updatedEvents)
          );
        } catch (e) {
          console.error(
            `[ResultTab:${sequenceIdx}] Failed to save SSE events to sessionStorage:`,
            e
          );
        }
        return updatedEvents;
      });
    },
    [sequenceIdx, jobId]
  );

  // --- Hook Usage ---
  useSSE(jobId, sequenceIdx, processSseEvent);

  // --- Mount/Unmount Effect ---
  useEffect(() => {
    isMounted.current = true;
    console.log(
      `[ResultTab:${sequenceIdx}] Component Did Mount. Initial state length: ${rawSseEvents.length}`
    );
    return () => {
      isMounted.current = false;
      console.log(`[ResultTab:${sequenceIdx}] Component Will Unmount.`);
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [sequenceIdx]);

  // --- Rendering ---
  return (
    <div className="sequence-results">
      {/* Status Messages */}
      {rawSseEvents.length === 0 &&
        !streamClosed.current &&
        !hasReceivedData.current && (
          <p className="text-gray-500 dark:text-gray-400 mb-4">
            Connecting to event stream...
          </p>
        )}

      {/* Protocol Tracker Component */}
      {protocolSteps.length > 0 && (
        <div className="border rounded-lg border-gray-200 dark:border-gray-700 overflow-hidden">
          <ProtocolTracker
            steps={protocolSteps}
            messages={protocolMessages}
            callouts={protocolCallouts}
            sseData={protocolSseData}
          />
        </div>
      )}
    </div>
  );
};

export default ResultTab;
