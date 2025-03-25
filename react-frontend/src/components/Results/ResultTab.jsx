import React, { useCallback, useState, useEffect } from "react";
import RestrictionSiteSummary from "./RestrictionSiteSummary";
import useSSE from "../../hooks/useSSE"; // adjust the path if needed

// Helper to format primer sequences consistently.
const formatPrimerSequence = (primer) => {
  if (!primer) return "None";
  if (typeof primer === "string") return primer;
  if (Array.isArray(primer)) return primer.join(", ");
  if (primer.sequence) return primer.sequence;
  return "None";
};

const ResultTab = ({ result, sequenceIdx }) => {
  const [copied, setCopied] = useState(false);
  // Initialize local progress state; default to 0% if not provided.
  const [progress, setProgress] = useState(
    result.progress || { percentage: 0, message: "" }
  );
  // New state to hold the restriction sites data.
  const [restrictionSites, setRestrictionSites] = useState(
    result.restriction_sites || []
  );

  const jobId = sessionStorage.getItem("jobId") || "";
  // Subscribe to SSE updates for this specific sequence.
  console.log(
    "Subscribing to SSE for jobId:",
    jobId,
    "sequenceIdx:",
    sequenceIdx
  );

  // Register this tab to receive tab-specific updates from the server using the jobId and sequenceIdx.
  const sseResult = useSSE(jobId, sequenceIdx);

  useEffect(() => {
    if (sseResult && sseResult.data) {
      const sseData = sseResult.data;
      setProgress({ percentage: sseData.progress, message: sseData.message });

      if (sseData.step === "Restriction Site Detection" && sseData.sites) {
        const sites = sseData.sites.map((site) => ({
          enzyme: site.enzyme,
          sequence: site.recognitionSeq,
          position: site.position,
          strand: site.strand,
        }));
        setRestrictionSites(sites);
      }
    }
  }, [sseResult]);

  // Copy PCR primer data to clipboard.
  const copyPrimersToClipboard = useCallback(() => {
    if (!result?.PCR_reactions) return;

    const rows = Object.entries(result.PCR_reactions).flatMap(
      ([reactionName, primers]) => {
        const forwardSeq = formatPrimerSequence(primers.forward);
        const reverseSeq = formatPrimerSequence(primers.reverse);
        return [
          `${reactionName}_FWD\t${forwardSeq}`,
          `${reactionName}_REV\t${reverseSeq}`,
        ];
      }
    );
    const finalText = rows.join("\n");

    navigator.clipboard
      .writeText(finalText)
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
      })
      .catch((err) => console.error("Failed to copy primers:", err));
  }, [result]);

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

  // Render PCR reactions if available.
  const renderPCRReactions = () => {
    if (result.PCR_reactions && Object.keys(result.PCR_reactions).length > 0) {
      return (
        <div className="pcr-summary section-container">
          <div className="section-header">
            <h3>PCR Reactions</h3>
            <button onClick={copyPrimersToClipboard} className="small-button">
              {copied ? "Copied!" : "Copy Primers"}
            </button>
          </div>
          <div className="table-container">
            <table>
              <thead>
                <tr>
                  <th>Reaction</th>
                  <th>Forward Primer</th>
                  <th>Reverse Primer</th>
                </tr>
              </thead>
              <tbody>
                {Object.entries(result.PCR_reactions).map(
                  ([reaction, primers], idx) => (
                    <tr key={idx}>
                      <td>{reaction}</td>
                      <td className="primer-cell">
                        {formatPrimerSequence(primers.forward)}
                      </td>
                      <td className="primer-cell">
                        {formatPrimerSequence(primers.reverse)}
                      </td>
                    </tr>
                  )
                )}
              </tbody>
            </table>
          </div>
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
      {renderPCRReactions()}
      {result.errors && (
        <div className="error-message">
          <strong>Error:</strong> {result.errors}
        </div>
      )}
    </div>
  );
};

export default ResultTab;
