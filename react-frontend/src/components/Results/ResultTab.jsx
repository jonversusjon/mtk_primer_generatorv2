import React, {
  useState,
  useEffect,
  useCallback,
  useRef,
  useMemo,
} from "react";
import useSSE from "../../hooks/useSSE";
import ProtocolTracker from "./ProtocolTracker";

const ResultTab = ({ result, sequenceIdx }) => {
  // Define all possible steps in order - ensure these match EXACTLY with backend step names
  const allSteps = [
    "Preprocessing",
    "Restriction Sites",
    "Mutation Analysis",
    "Primer Design",
    "PCR Reaction Grouping",
  ];

  // Initialize steps with waiting status
  const [protocolSteps, setProtocolSteps] = useState(
    allSteps.map((name) => ({
      name,
      status: "waiting",
      progress: null,
      message: "",
      notificationCount: 0,
    }))
  );

  // Use a ref to track processed events to avoid duplicates
  const processedEvents = useRef(new Set());

  // Messages for the message log - use callback for initialization
  const [messagesSet, setMessagesSet] = useState(() => new Set());

  // Create a map to convert step names to data keys
  const stepToDataKeyMap = useMemo(
    () => ({
      Preprocessing: "Preprocessing",
      "Restriction Site Detection": "RestrictionSiteDetection",
      "Mutation Analysis": "MutationAnalysis",
      "Primer Design": "PrimerDesign",
      "PCR Reaction Grouping": "PCRReactionGrouping",
    }),
    []
  );

  // State for result data organized by step
  const [stepData, setStepData] = useState({
    Preprocessing: {
      processedSequence: result.processed_sequence || "",
    },
    RestrictionSiteDetection: {
      restrictionSites: result.restriction_sites || [],
    },
    MutationAnalysis: {
      mutations: [],
    },
    PrimerDesign: {
      edgePrimers: result.edge_primers || null,
      mutPrimers: result.mut_primers || {},
    },
    PCRReactionGrouping: {
      pcrReactions: result.PCR_reactions || [],
    },
  });

  // Store SSE data per step instead of just the latest overall
  const [sseDataByStep, setSseDataByStep] = useState({});

  const jobId = sessionStorage.getItem("jobId") || "";

  // Register this tab to receive tab-specific updates
  const sseResult = useSSE(jobId, sequenceIdx);

  useEffect(() => {
    return () =>
      console.log(
        "Component using useSSE for",
        jobId,
        sequenceIdx,
        "is unmounting"
      );
  }, [jobId, sequenceIdx]);

  // Debug log all SSE events
  useEffect(() => {
    if (sseResult) {
      console.log(`[ResultTab:${sequenceIdx}] SSE event received:`, sseResult);
    }
  }, [sseResult, sequenceIdx]);

  // Update step data based on SSE event
  const updateStepData = useCallback(
    (sseData) => {
      // Skip if step is missing
      if (!sseData.step) {
        console.warn(
          `[ResultTab:${sequenceIdx}] Missing step in SSE data:`,
          sseData
        );
        return;
      }

      // Get the data key for this step
      const dataKey = stepToDataKeyMap[sseData.step];
      if (!dataKey) {
        console.warn(
          `[ResultTab:${sequenceIdx}] No data key mapping for step:`,
          sseData.step
        );
        return;
      }

      switch (sseData.step) {
        case "Restriction Sites":
          if (sseData.sites) {
            const sites = sseData.sites.map((site) => ({
              enzyme: site.enzyme,
              recognition_seq: site.recognitionSeq,
              position: site.position,
              strand: site.strand,
            }));

            setStepData((prevData) => ({
              ...prevData,
              RestrictionSiteDetection: {
                ...prevData.RestrictionSiteDetection,
                restrictionSites: sites,
              },
            }));

            // Update notification count for this step
            if (sites.length > 0) {
              setProtocolSteps((prevSteps) => {
                return prevSteps.map((step) =>
                  step.name === "Restriction Sites"
                    ? { ...step, notificationCount: sites.length }
                    : step
                );
              });
            }
          }
          break;

        case "Mutation Analysis":
          if (sseData.mutations) {
            const mutations = sseData.mutations.map((mutation) => ({
              type: mutation.type,
              position: mutation.position,
              sequence: mutation.sequence,
            }));

            setStepData((prevData) => ({
              ...prevData,
              MutationAnalysis: {
                ...prevData.MutationAnalysis,
                mutations: mutations,
              },
            }));

            // Update notification count for this step
            if (mutations.length > 0) {
              setProtocolSteps((prevSteps) => {
                return prevSteps.map((step) =>
                  step.name === "Mutation Analysis"
                    ? { ...step, notificationCount: mutations.length }
                    : step
                );
              });
            }
          }
          break;
        case "Primer Design":
          // Track primer updates for notification count
          let edgePrimerCount = 0;
          let mutPrimerCount = 0;

          if (sseData.edgePrimers) {
            edgePrimerCount = Object.keys(sseData.edgePrimers).length;
            setStepData((prevData) => ({
              ...prevData,
              PrimerDesign: {
                ...prevData.PrimerDesign,
                edgePrimers: sseData.edgePrimers,
              },
            }));
          }

          if (sseData.mutPrimers) {
            mutPrimerCount = Object.keys(sseData.mutPrimers).length;
            setStepData((prevData) => ({
              ...prevData,
              PrimerDesign: {
                ...prevData.PrimerDesign,
                mutPrimers: sseData.mutPrimers,
              },
            }));
          }

          // Update notification count if we have any primers
          const totalPrimerCount = edgePrimerCount + mutPrimerCount;
          if (totalPrimerCount > 0) {
            setProtocolSteps((prevSteps) => {
              return prevSteps.map((step) =>
                step.name === "Primer Design"
                  ? { ...step, notificationCount: totalPrimerCount }
                  : step
              );
            });
          }
          break;

        case "PCR Reaction Grouping":
          if (
            sseData.domestication_result &&
            sseData.domestication_result.pcr_reactions
          ) {
            const pcrReactions = sseData.domestication_result.pcr_reactions;
            setStepData((prevData) => ({
              ...prevData,
              PCRReactionGrouping: {
                ...prevData.PCRReactionGrouping,
                pcrReactions: pcrReactions,
              },
            }));

            // Update notification count for PCR reactions
            if (pcrReactions.length > 0) {
              setProtocolSteps((prevSteps) => {
                return prevSteps.map((step) =>
                  step.name === "PCR Reaction Grouping"
                    ? { ...step, notificationCount: pcrReactions.length }
                    : step
                );
              });
            }
          }
          break;

        case "Preprocessing":
          console.log(
            `[ResultTab:${sequenceIdx}] Preprocessing step received:`,
            sseData
          );
          if (sseData.processedSequence) {
            console.log(
              `[ResultTab:${sequenceIdx}] Updating processed sequence:`,
              sseData.processedSequence
            );
            setStepData((prevData) => ({
              ...prevData,
              Preprocessing: {
                ...prevData.Preprocessing,
                processedSequence: sseData.processedSequence,
              },
            }));
          } else {
            console.warn(
              `[ResultTab:${sequenceIdx}] No processedSequence found in SSE data.`
            );
          }
          break;

        default:
          console.warn(
            `[ResultTab:${sequenceIdx}] Unknown step:`,
            sseData.step
          );
          break;
      }
    },
    [sequenceIdx, stepToDataKeyMap]
  );

  // Process SSE data - memoized to ensure consistent reference
  const processSseData = useCallback(
    (sseData) => {
      if (!sseData || !sseData.step) return;

      // Generate a unique ID for this event - include timestamp if available
      const eventId = `${sseData.step}-${sseData.message}-${
        sseData.stepProgress
      }-${Date.now()}`;
      if (processedEvents.current.has(eventId)) {
        console.log(
          `[ResultTab:${sequenceIdx}] Skipping duplicate event:`,
          eventId
        );
        return;
      }

      // Mark as processed
      processedEvents.current.add(eventId);
      console.log(
        `[ResultTab:${sequenceIdx}] Processing event:`,
        eventId,
        sseData
      );

      // Update SSE data by step - store each step's data separately
      setSseDataByStep((prevData) => ({
        ...prevData,
        [sseData.step]: sseData,
      }));

      // Update steps state
      setProtocolSteps((prevSteps) => {
        const stepIndex = prevSteps.findIndex(
          (step) => step.name === sseData.step
        );

        if (stepIndex === -1) {
          console.warn(
            `[ResultTab:${sequenceIdx}] Step not found:`,
            sseData.step
          );
          return prevSteps; // Step not found, return unchanged
        }

        // Create a new steps array to modify
        const newSteps = [...prevSteps];

        // Mark all previous steps as completed if they're not already
        for (let i = 0; i < stepIndex; i++) {
          if (newSteps[i].status !== "completed") {
            newSteps[i] = {
              ...newSteps[i],
              status: "completed",
              progress: 100,
            };
          }
        }

        // Update the current step
        const stepProgress =
          sseData.stepProgress !== undefined
            ? sseData.stepProgress
            : newSteps[stepIndex].progress;
        const stepMessage = sseData.message || newSteps[stepIndex].message;

        if (stepProgress >= 100) {
          // If step reached 100%, mark as completed
          newSteps[stepIndex] = {
            ...newSteps[stepIndex],
            status: "completed",
            progress: 100,
            message: stepMessage,
          };

          // If this wasn't the last step, set the next step to active
          if (stepIndex < newSteps.length - 1) {
            newSteps[stepIndex + 1] = {
              ...newSteps[stepIndex + 1],
              status: "active",
            };
          }
        } else {
          // This is an in-progress update
          newSteps[stepIndex] = {
            ...newSteps[stepIndex],
            status: "active",
            progress: stepProgress,
            message: stepMessage,
          };
        }

        return newSteps;
      });

      // Add to messages log if there's a new message
      if (sseData.message) {
        setMessagesSet((prev) => {
          // Create a new Set to trigger a re-render
          const newSet = new Set(prev);
          newSet.add(`${sseData.step}: ${sseData.message}`);
          return newSet;
        });
      }

      // Update step-specific data based on event content
      updateStepData(sseData);
    },
    [sequenceIdx, updateStepData]
  );

  // Handle SSE updates - now with proper dependencies
  useEffect(() => {
    if (!sseResult) return;

    // Handle different data structures that might come from SSE
    let sseData;
    if (sseResult.data) {
      sseData = sseResult.data;
    } else {
      sseData = sseResult;
    }

    if (sseData) {
      console.log(`[ResultTab:${sequenceIdx}] Processing SSE Data:`, sseData);
      processSseData(sseData);
    }
  }, [sseResult, processSseData, sequenceIdx]);

  // Convert message set to array for rendering
  const messages = Array.from(messagesSet);

  return (
    <div className="sequence-results p-4">
      <div className="mb-6">
        <ProtocolTracker
          steps={protocolSteps}
          messages={messages}
          resultData={stepData}
          sseData={sseDataByStep}
        />
      </div>
    </div>
  );
};

export default ResultTab;
