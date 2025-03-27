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
      "Restriction Sites": "RestrictionSiteDetection",
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
      restrictionSites: result.restriction_sites || [], // Add this line to access restriction sites
      mutationSets: result.mut_primers
        ? [{ mutations: result.mut_primers }]
        : [], // Add this line for mutation sets
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

  // Helper functions for deep merging objects and handling arrays
  const mergeDeep = useCallback((target, source) => {
    if (!source) return target;

    const output = { ...target };

    Object.keys(source).forEach((key) => {
      if (source[key] === null || source[key] === undefined) {
        // Skip null/undefined values to preserve existing data
        return;
      }

      // Handle explicit empty arrays (clear the data)
      if (Array.isArray(source[key]) && source[key].length === 0) {
        output[key] = [];
        return;
      }

      // If both are objects and not arrays, recursively merge
      if (
        typeof source[key] === "object" &&
        source[key] !== null &&
        typeof output[key] === "object" &&
        output[key] !== null &&
        !Array.isArray(source[key]) &&
        !Array.isArray(output[key])
      ) {
        output[key] = mergeDeep(output[key], source[key]);
      }
      // For arrays, append new items instead of replacing (unless specified as empty)
      else if (Array.isArray(source[key]) && Array.isArray(output[key])) {
        // Check if arrays contain objects with IDs for deduplication
        if (
          source[key].length > 0 &&
          typeof source[key][0] === "object" &&
          source[key][0] !== null
        ) {
          // If objects have an ID field, use that for deduplication
          const idField =
            "id" in source[key][0]
              ? "id"
              : "enzyme" in source[key][0]
              ? "enzyme"
              : "position" in source[key][0]
              ? "position"
              : null;

          if (idField) {
            // Filter out existing items with the same ID
            const existingIds = new Set(
              output[key].map((item) => item[idField])
            );
            const newItems = source[key].filter(
              (item) => !existingIds.has(item[idField])
            );
            output[key] = [...output[key], ...newItems];
          } else {
            // No ID field for deduplication, just append
            output[key] = [...output[key], ...source[key]];
          }
        } else {
          // Simple values, just append
          output[key] = [...output[key], ...source[key]];
        }
      }
      // Otherwise just replace the value
      else {
        output[key] = source[key];
      }
    });

    return output;
  }, []);

  // Update step data based on SSE event - modified for incremental updates
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

            setStepData((prevData) => {
              // Get existing sites or empty array
              const existingSites =
                prevData.RestrictionSiteDetection?.restrictionSites || [];

              // Deduplicate by enzyme and position
              const existingKeys = new Set(
                existingSites.map(
                  (site) => `${site.enzyme}-${site.position}-${site.strand}`
                )
              );

              // Filter out duplicates
              const newSites = sites.filter(
                (site) =>
                  !existingKeys.has(
                    `${site.enzyme}-${site.position}-${site.strand}`
                  )
              );

              // Create combined array
              const combinedSites = [...existingSites, ...newSites];

              // Update notification count for combined sites
              if (combinedSites.length > 0) {
                setProtocolSteps((prevSteps) => {
                  return prevSteps.map((step) =>
                    step.name === "Restriction Sites"
                      ? { ...step, notificationCount: combinedSites.length }
                      : step
                  );
                });
              }

              return {
                ...prevData,
                RestrictionSiteDetection: {
                  ...prevData.RestrictionSiteDetection,
                  restrictionSites: combinedSites,
                },
              };
            });
          }
          break;

        case "Mutation Analysis":
          if (sseData.mutations) {
            const mutations = sseData.mutations.map((mutation) => ({
              type: mutation.type,
              position: mutation.position,
              sequence: mutation.sequence,
            }));

            setStepData((prevData) => {
              // Get existing mutations or empty array
              const existingMutations =
                prevData.MutationAnalysis?.mutations || [];

              // Deduplicate by position and type
              const existingKeys = new Set(
                existingMutations.map((mut) => `${mut.type}-${mut.position}`)
              );

              // Filter out duplicates
              const newMutations = mutations.filter(
                (mut) => !existingKeys.has(`${mut.type}-${mut.position}`)
              );

              // Create combined array
              const combinedMutations = [...existingMutations, ...newMutations];

              // Update notification count for combined mutations
              if (combinedMutations.length > 0) {
                setProtocolSteps((prevSteps) => {
                  return prevSteps.map((step) =>
                    step.name === "Mutation Analysis"
                      ? { ...step, notificationCount: combinedMutations.length }
                      : step
                  );
                });
              }

              // Handle mutation sets if provided in the SSE data
              let updatedMutationSets =
                prevData.MutationAnalysis?.mutationSets || [];
              if (sseData.mutationSets) {
                updatedMutationSets = sseData.mutationSets;
              } else if (sseData.mut_primers) {
                // Convert mut_primers to mutation sets format if needed
                updatedMutationSets = [{ mutations: sseData.mut_primers }];
              }

              // Get restriction site data if provided or use existing
              const restrictionSites =
                sseData.restrictionSites ||
                prevData.MutationAnalysis?.restrictionSites ||
                prevData.RestrictionSiteDetection?.restrictionSites ||
                [];

              return {
                ...prevData,
                MutationAnalysis: {
                  ...prevData.MutationAnalysis,
                  mutations: combinedMutations,
                  restrictionSites: restrictionSites,
                  mutationSets: updatedMutationSets,
                },
              };
            });
          }
          break;

        case "Primer Design":
          setStepData((prevData) => {
            // Starting with the existing data
            const updatedPrimerData = {
              ...prevData.PrimerDesign,
            };

            // Update edgePrimers if provided
            if (sseData.edgePrimers) {
              // Merge with existing edge primers rather than replacing
              updatedPrimerData.edgePrimers = {
                ...(updatedPrimerData.edgePrimers || {}),
                ...sseData.edgePrimers,
              };
            }

            // Update mutPrimers if provided
            if (sseData.mutPrimers) {
              // Merge with existing mutation primers rather than replacing
              updatedPrimerData.mutPrimers = {
                ...(updatedPrimerData.mutPrimers || {}),
                ...sseData.mutPrimers,
              };
            }

            // Calculate total primer count for notification
            const edgePrimerCount = updatedPrimerData.edgePrimers
              ? Object.keys(updatedPrimerData.edgePrimers).length
              : 0;

            const mutPrimerCount = updatedPrimerData.mutPrimers
              ? Object.keys(updatedPrimerData.mutPrimers).length
              : 0;

            const totalPrimerCount = edgePrimerCount + mutPrimerCount;

            // Update notification count if we have any primers
            if (totalPrimerCount > 0) {
              setProtocolSteps((prevSteps) => {
                return prevSteps.map((step) =>
                  step.name === "Primer Design"
                    ? { ...step, notificationCount: totalPrimerCount }
                    : step
                );
              });
            }

            return {
              ...prevData,
              PrimerDesign: updatedPrimerData,
            };
          });
          break;

        case "PCR Reaction Grouping":
          if (
            sseData.domestication_result &&
            sseData.domestication_result.pcr_reactions
          ) {
            const newPcrReactions = sseData.domestication_result.pcr_reactions;

            setStepData((prevData) => {
              // Get existing PCR reactions or empty array
              const existingReactions =
                prevData.PCRReactionGrouping?.pcrReactions || [];

              // For PCR reactions, we need to check if the reaction is already present
              // We'll use a combination of template, forwardPrimer and reversePrimer as a key
              const existingKeys = new Set(
                existingReactions.map(
                  (rxn) =>
                    `${rxn.template || ""}-${rxn.forwardPrimer || ""}-${
                      rxn.reversePrimer || ""
                    }`
                )
              );

              // Filter out duplicates
              const uniqueNewReactions = newPcrReactions.filter(
                (rxn) =>
                  !existingKeys.has(
                    `${rxn.template || ""}-${rxn.forwardPrimer || ""}-${
                      rxn.reversePrimer || ""
                    }`
                  )
              );

              // Create combined array
              const combinedReactions = [
                ...existingReactions,
                ...uniqueNewReactions,
              ];

              // Update notification count
              if (combinedReactions.length > 0) {
                setProtocolSteps((prevSteps) => {
                  return prevSteps.map((step) =>
                    step.name === "PCR Reaction Grouping"
                      ? { ...step, notificationCount: combinedReactions.length }
                      : step
                  );
                });
              }

              return {
                ...prevData,
                PCRReactionGrouping: {
                  ...prevData.PCRReactionGrouping,
                  pcrReactions: combinedReactions,
                },
              };
            });
          }
          break;

          case "Preprocessing":
            console.log(
              `[ResultTab:${sequenceIdx}] Preprocessing step received:`,
              sseData
            );
          
            setStepData((prevData) => ({
              ...prevData,
              Preprocessing: {
                ...prevData.Preprocessing,
                ...(sseData.processedSequence && {
                  processedSequence: sseData.processedSequence,
                }),
              },
            }));
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

      if (sseData.sequenceIdx !== sequenceIdx) {
        console.log(`[ResultTab:${sequenceIdx}] Ignoring event for sequenceIdx ${sseData.sequenceIdx}`);
        return;
      }
      // Generate a unique ID for this event - include timestamp if available
      const eventId = `${sseData.sequenceIdx}-${sseData.step}-${sseData.message}-${sseData.stepProgress}-${Date.now()}`;
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
      // Using mergeDeep to combine with existing data instead of replacing
      setSseDataByStep((prevData) => {
        // If this step already has data, merge with it
        if (prevData[sseData.step]) {
          return {
            ...prevData,
            [sseData.step]: mergeDeep(prevData[sseData.step], sseData),
          };
        }
        // Otherwise just add the new data
        return {
          ...prevData,
          [sseData.step]: sseData,
        };
      });
  
      // 3) If SSE includes a notification_count, apply it
      if (sseData.notification_count && sseData.notification_count > 0) {
        setProtocolSteps((prevSteps) =>
          prevSteps.map((step) =>
            step.name === sseData.step
              ? { ...step, notificationCount: sseData.notification_count }
              : step
          )
        );
      }
  
      // 4) Update step status and data (progress, message, arrays, etc.)
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
    [sequenceIdx, updateStepData, mergeDeep]
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
