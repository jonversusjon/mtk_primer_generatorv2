import React, { useState, useEffect, useCallback, useRef } from "react";
import useSSE from "../../hooks/useSSE";
import ProtocolTracker from "./ProtocolTracker";
import {mergeDeep} from "../../utils/mergeDeep";

const ResultTab = ({ jobId, sequenceIdx }) => {
  // This will hold all SSE data keyed by step name:
  // e.g. sseDataByStep["Restriction Sites"] = { sites: [...], ... }
  const [sseDataByStep, setSseDataByStep] = useState({});

  // The steps array for the ProtocolTracker
  const [protocolSteps, setProtocolSteps] = useState([
    {
      name: "Preprocessing",
      status: "waiting",
      progress: 0,
      message: "",
      notificationCount: 0,
    },
    {
      name: "Restriction Sites",
      status: "waiting",
      progress: 0,
      message: "",
      notificationCount: 0,
    },
    {
      name: "Mutation Analysis",
      status: "waiting",
      progress: 0,
      message: "",
      notificationCount: 0,
    },
    {
      name: "Primer Design",
      status: "waiting",
      progress: 0,
      message: "",
      notificationCount: 0,
    },
    {
      name: "PCR Reaction Grouping",
      status: "waiting",
      progress: 0,
      message: "",
      notificationCount: 0,
    },
  ]);

  // A set for deduplicating SSE events
  const processedEvents = useRef(new Set());

  // A set for collecting text messages (to show in the “messages” log)
  const [messagesSet, setMessagesSet] = useState(new Set());

  // Connect to SSE (custom hook)
  const sseResult = useSSE(jobId, sequenceIdx);

  // Process SSE events
  const processSseData = useCallback(
    (sseData) => {
      if (!sseData || !sseData.step) return;

      // Filter out events for other sequences
      if (sseData.sequenceIdx !== sequenceIdx) {
        console.log(
          `[ResultTab:${sequenceIdx}] Ignoring event for sequenceIdx ${sseData.sequenceIdx}`
        );

        return;
      }

      // Deduplicate events
      const eventId = `${sseData.sequenceIdx}-${sseData.step}-${sseData.message}-${sseData.stepProgress}`;
      if (processedEvents.current.has(eventId)) return;
      processedEvents.current.add(eventId);

      // Merge partial SSE data into our single sseDataByStep object
      setSseDataByStep((prev) => {
        const currentForStep = prev[sseData.step] || {};
        const merged = mergeDeep(currentForStep, sseData);
        return { ...prev, [sseData.step]: merged };
      });

      // If SSE includes a notification_count, update protocolSteps
      if (sseData.notification_count > 0) {
        setProtocolSteps((prevSteps) =>
          prevSteps.map((step) =>
            step.name === sseData.step
              ? { ...step, notificationCount: sseData.notification_count }
              : step
          )
        );
      }

      // Update step status, progress, etc.
      setProtocolSteps((prevSteps) => {
        const stepIndex = prevSteps.findIndex((s) => s.name === sseData.step);
        if (stepIndex < 0) return prevSteps; // not found
        const newSteps = [...prevSteps];

        // Mark previous steps completed
        for (let i = 0; i < stepIndex; i++) {
          if (newSteps[i].status !== "completed") {
            newSteps[i] = {
              ...newSteps[i],
              status: "completed",
              progress: 100,
            };
          }
        }

        // This step’s progress
        const stepProgress =
          sseData.stepProgress ?? newSteps[stepIndex].progress;
        const stepMessage = sseData.message ?? newSteps[stepIndex].message;

        if (stepProgress >= 100) {
          newSteps[stepIndex] = {
            ...newSteps[stepIndex],
            status: "completed",
            progress: 100,
            message: stepMessage,
          };
          // Make the next step active if it exists
          if (stepIndex < newSteps.length - 1) {
            newSteps[stepIndex + 1] = {
              ...newSteps[stepIndex + 1],
              status: "active",
            };
          }
        } else {
          // Mark this step as active/in-progress
          newSteps[stepIndex] = {
            ...newSteps[stepIndex],
            status: "active",
            progress: stepProgress,
            message: stepMessage,
          };
        }

        return newSteps;
      });

      // Add text message to messagesSet
      if (sseData.message) {
        setMessagesSet((prev) => {
          const newSet = new Set(prev);
          newSet.add(`${sseData.step}: ${sseData.message}`);
          return newSet;
        });
      }
    },
    [sequenceIdx]
  );

  // Listen for SSE events
  useEffect(() => {
    if (sseResult?.data) {
      processSseData(sseResult.data);
    }
    // Or if sseResult is the raw SSE object itself:
    else if (sseResult?.step) {
      processSseData(sseResult);
    }
  }, [sseResult, processSseData]);

  // Convert messagesSet to array
  const messages = Array.from(messagesSet);

  return (
    <div className="sequence-results p-4">
      <ProtocolTracker
        steps={protocolSteps}
        messages={messages}
        sseData={sseDataByStep} // Pass all SSE data
      />
    </div>
  );
};

export default ResultTab;
