import React, { useState, useEffect } from "react";
import RestrictionSiteSummary from "./RestrictionSiteSummary";
import useSSE from "../../hooks/useSSE";
import PcrReaction from "./PcrReactions"; // Adjust the path as needed

const ResultTab = ({ result, sequenceIdx }) => {
  // Initialize local progress state; default to 0% if not provided.
  const [progress, setProgress] = useState(
    result.progress || { percentage: 0, message: "" }
  );

  // New state to hold the restriction sites data.
  const [restrictionSites, setRestrictionSites] = useState(
    result.restriction_sites || []
  );
  const [reactions, setReactions] = useState(result.reactions || []);

  const jobId = sessionStorage.getItem("jobId") || "";

  // Register this tab to receive tab-specific updates from the server using the jobId and sequenceIdx.
  const sseResult = useSSE(jobId, sequenceIdx);

  useEffect(() => {
    if (sseResult) {
      console.log("SSE Data Received:", sseResult);
      if (sseResult.data) {
        const sseData = sseResult.data;

        // Update progress state regardless of event type
        if (sseData.progress !== undefined && sseData.message) {
          setProgress({
            percentage: sseData.progress,
            message: sseData.message,
          });
        }

        // Handle events based on the step, regardless of type
        switch (sseData.step) {
          case "Restriction Site Detection":
            if (sseData.sites) {
              const sites = sseData.sites.map((site) => ({
                enzyme: site.enzyme,
                recognition_seq: site.recognitionSeq,
                position: site.position,
                strand: site.strand,
              }));
              setRestrictionSites(sites);
            }
            break;

          case "PCR Reaction Grouping":
            // Check both for data or progress events
            if (
              sseData.domestication_result &&
              sseData.domestication_result.pcr_reactions
            ) {
              setReactions(sseData.domestication_result.pcr_reactions);
              console.log(
                "Updated reactions:",
                sseData.domestication_result.pcr_reactions
              );
            }
            break;

          // Add additional cases as needed.
          default:
            break;
        }
      }
    }
  }, [sseResult]);

  // Render a progress bar if the process isn’t complete.
  const renderProgress = () => {
    if (progress && progress.percentage < 100) {
      return (
        <div className="progress-container">
          <div
            className="progress-bar"
            style={{ width: `${progress.percentage}%` }}
          ></div>
          <div className="progress-message">{progress.message}</div>
        </div>
      );
    }
    return null;
  };

  // Render a placeholder if the sequence is still processing.
  const renderPlaceholderMessage = () => {
    if (
      result.placeholder &&
      !result.PCR_reactions &&
      (!progress || progress.percentage < 100)
    ) {
      return (
        <div className="placeholder-message">
          This sequence is still processing...
        </div>
      );
    }
    return null;
  };

  return (
    <div className="sequence-results">
      {renderProgress()}
      {renderPlaceholderMessage()}
      {result.messages && result.messages.length > 0 && (
        <div className="messages">
          {result.messages.map((msg, idx) => (
            <div key={idx}>{msg}</div>
          ))}
        </div>
      )}
      <div className="mtk-part-info">
        {result.mtk_part_left === result.mtk_part_right ? (
          <p>
            <strong>MTK Part Number:</strong> {result.mtk_part_left}
          </p>
        ) : (
          <p>
            <strong>MTK Part Number Left:</strong> {result.mtk_part_left} <br />
            <strong>MTK Part Number Right:</strong> {result.mtk_part_right}
          </p>
        )}
      </div>
      {/* Render the updated restriction sites */}
      {restrictionSites && restrictionSites.length > 0 && (
        <RestrictionSiteSummary sites={restrictionSites} />
      )}
      {/* Render PCR reactions via the new component */}
      <PcrReaction pcrReactions={reactions} />
      {result.errors && (
        <div className="error-message">
          <strong>Error:</strong> {result.errors}
        </div>
      )}
    </div>
  );
};

export default ResultTab;
