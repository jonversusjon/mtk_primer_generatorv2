import React, { useState, useEffect, useCallback, useRef } from "react";
import useSSE from "../../hooks/useSSE";
import ProtocolTracker from "./ProtocolTracker";

const ResultTab = ({ jobId, sequenceIdx }) => {
  // This state holds SSE data for each step keyed by step name.
  const [sseDataByStep, setSseDataByStep] = useState({});

  // Local model for protocol steps – these drive progress and notifications.
  const [protocolSteps, setProtocolSteps] = useState([
    { name: "Preprocessing", status: "waiting", progress: 0, message: "", notificationCount: 0 },
    { name: "Restriction Sites", status: "waiting", progress: 0, message: "", notificationCount: 0 },
    { name: "Mutation Analysis", status: "waiting", progress: 0, message: "", notificationCount: 0 },
    { name: "Primer Design", status: "waiting", progress: 0, message: "", notificationCount: 0 },
    { name: "PCR Reaction Grouping", status: "waiting", progress: 0, message: "", notificationCount: 0 },
  ]);

  // A Set to deduplicate SSE events.
  const processedEvents = useRef(new Set());
  // A Set to collect messages for callouts.
  const [messagesSet, setMessagesSet] = useState(new Set());

  // Subscribe to SSE using our custom hook.
  const sseResult = useSSE(jobId, sequenceIdx);

  const processSseData = useCallback(
    (sseData) => {
      if (!sseData) return;

      // If this is a final result event (detected by presence of "result" property),
      // log it and ignore it (do not update state).
      if (sseData.result !== undefined) {
        console.log(`[ResultTab:${sequenceIdx}] Final result received:`, sseData);
        return;
      }

      // Ensure the event has a "step" property.
      if (!sseData.step) return;

      // Only process events for the correct sequence.
      if (sseData.sequenceIdx !== sequenceIdx) {
        console.log(`[ResultTab:${sequenceIdx}] Ignoring event for sequenceIdx ${sseData.sequenceIdx}`);
        return;
      }

      // Generate a unique event id to deduplicate.
      const eventId = `${sseData.sequenceIdx}-${sseData.step}-${sseData.message}-${sseData.stepProgress}`;
      if (processedEvents.current.has(eventId)) return;
      processedEvents.current.add(eventId);

      // Update the SSE data for the step.
      setSseDataByStep((prev) => ({
        ...prev,
        [sseData.step]: sseData,
      }));

      // Update protocol steps if a notification count is provided.
      if (sseData.notification_count > 0) {
        setProtocolSteps((prevSteps) =>
          prevSteps.map((step) =>
            step.name === sseData.step
              ? { ...step, notificationCount: sseData.notification_count }
              : step
          )
        );
      }

      // Update progress and status for the step.
      setProtocolSteps((prevSteps) => {
        const stepIndex = prevSteps.findIndex((s) => s.name === sseData.step);
        if (stepIndex < 0) return prevSteps;
        const newSteps = [...prevSteps];

        // Mark all previous steps as completed.
        for (let i = 0; i < stepIndex; i++) {
          if (newSteps[i].status !== "completed") {
            newSteps[i] = { ...newSteps[i], status: "completed", progress: 100 };
          }
        }

        const stepProgress = sseData.stepProgress ?? newSteps[stepIndex].progress;
        const stepMessage = sseData.message ?? newSteps[stepIndex].message;

        if (stepProgress >= 100) {
          newSteps[stepIndex] = { ...newSteps[stepIndex], status: "completed", progress: 100, message: stepMessage };
          if (stepIndex < newSteps.length - 1) {
            newSteps[stepIndex + 1] = { ...newSteps[stepIndex + 1], status: "active" };
          }
        } else {
          newSteps[stepIndex] = { ...newSteps[stepIndex], status: "active", progress: stepProgress, message: stepMessage };
        }
        return newSteps;
      });

      // Record the message in the callouts log.
      if (sseData.message) {
        setMessagesSet((prev) => {
          const newSet = new Set(prev);
          newSet.add({ step: sseData.step, message: sseData.message });
          return newSet;
        });
      }
    },
    [sequenceIdx]
  );

  useEffect(() => {
    if (sseResult) {
      if (sseResult.data) {
        console.log(`[ResultTab:${sequenceIdx}] Received SSE result (data):`, sseResult.data);
        processSseData(sseResult.data);
      } else {
        console.log(`[ResultTab:${sequenceIdx}] Received SSE result (raw):`, sseResult);
        processSseData(sseResult);
      }
    }
  }, [sseResult, processSseData, sequenceIdx]);

  // Convert messagesSet to an array.
  const messages = Array.from(messagesSet);

  return (
    <div className="sequence-results p-4">
      <ProtocolTracker steps={protocolSteps} messages={messages} sseData={sseDataByStep} />
    </div>
  );
};

export default ResultTab;
