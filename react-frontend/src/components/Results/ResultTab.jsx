import React, { useCallback, useState } from "react";
import RestrictionSiteSummary from "./RestrictionSiteSummary";

/*
  ResultTab displays the detailed information for a single sequence.
  - It renders a progress bar if the sequence is still processing.
  - It shows a placeholder message if the sequence is marked as a placeholder and has no final data.
  - It displays final results (e.g., PCR reactions, MTK parts, messages) when available.
*/

// Helper to format primer sequences consistently.
const formatPrimerSequence = (primer) => {
  if (!primer) return "None";
  if (typeof primer === "string") return primer;
  if (Array.isArray(primer)) return primer.join(", ");
  if (primer.sequence) return primer.sequence;
  return "None";
};

const ResultTab = ({ result }) => {
  const [copied, setCopied] = useState(false);
  const { mtk_part_left = "N/A", mtk_part_right = "N/A" } = result || {};

  // Copy PCR primer data (tab-separated) to the clipboard.
  const copyPrimersToClipboard = useCallback(() => {
    if (!result?.PCR_reactions) return;

    // Build rows for each PCR reaction.
    const rows = Object.entries(result.PCR_reactions).flatMap(([reactionName, primers]) => {
      const forwardSeq = formatPrimerSequence(primers.forward);
      const reverseSeq = formatPrimerSequence(primers.reverse);
      return [`${reactionName}_FWD\t${forwardSeq}`, `${reactionName}_REV\t${reverseSeq}`];
    });
    const finalText = rows.join("\n");

    navigator.clipboard
      .writeText(finalText)
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
      })
      .catch((err) => console.error("Failed to copy primers:", err));
  }, [result]);

  // Render progress bar if the sequence is still processing.
  const renderProgress = () => {
    if (result.progress && result.progress.percentage < 100) {
      return (
        <div className="progress-container">
          <div className="progress-bar" style={{ width: `${result.progress.percentage}%` }}></div>
          <div className="progress-message">{result.progress.message}</div>
        </div>
      );
    }
    return null;
  };

  // Render a placeholder message for sequences still processing and with no final data.
  const renderPlaceholderMessage = () => {
    if (result.placeholder && (!result.PCR_reactions && (!result.progress || result.progress.percentage < 100))) {
      return <div className="placeholder-message">This sequence is still processing...</div>;
    }
    return null;
  };

  // Render PCR reactions table if available.
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
                {Object.entries(result.PCR_reactions).map(([reaction, primers], idx) => (
                  <tr key={idx}>
                    <td>{reaction}</td>
                    <td className="primer-cell">{formatPrimerSequence(primers.forward)}</td>
                    <td className="primer-cell">{formatPrimerSequence(primers.reverse)}</td>
                  </tr>
                ))}
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
        {mtk_part_left === mtk_part_right ? (
          <p>
            <strong>MTK Part Number:</strong> {mtk_part_left}
          </p>
        ) : (
          <p>
            <strong>MTK Part Number Left:</strong> {mtk_part_left} <br />
            <strong>MTK Part Number Right:</strong> {mtk_part_right}
          </p>
        )}
      </div>
      {result.restriction_sites && result.restriction_sites.length > 0 && (
        <RestrictionSiteSummary sites={result.restriction_sites} />
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
