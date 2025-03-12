import React from "react";
import RestrictionSiteSummary from "./RestrictionSiteSummary";

function formatPrimerSequence(primer) {
  if (!primer) return "None";
  if (typeof primer === "string") return primer;
  if (Array.isArray(primer)) return primer.join(", ");
  if (primer.sequence) return primer.sequence;
  return "None";
}

function ResultTab({ result, index }) {
  const [copied, setCopied] = React.useState(false);

  const mtkPartLeft = result?.mtk_part_left || "N/A";
  const mtkPartRight = result?.mtk_part_right || "N/A";

  console.log("MTK Part Left:", mtkPartLeft);
  console.log("MTK Part Right:", mtkPartRight);

  // Copy all primer data (tab-separated) to the clipboard
  const copyPrimersToClipboard = React.useCallback(() => {
    if (!result || !result.PCR_reactions) {
      return;
    }

    const rows = [];
    // Iterate over each PCR reaction
    for (const [reactionName, primers] of Object.entries(
      result.PCR_reactions
    )) {
      // Forward primer
      const forwardSeq = formatPrimerSequence(primers.forward);
      rows.push(`${reactionName}_FWD\t${forwardSeq}\t`);
      // Reverse primer
      const reverseSeq = formatPrimerSequence(primers.reverse);
      rows.push(`${reactionName}_REV\t${reverseSeq}\t`);
    }

    const finalText = rows.join("\n");

    navigator.clipboard
      .writeText(finalText)
      .then(() => {
        setCopied(true);
        // Reset the button after 2 seconds
        setTimeout(() => {
          setCopied(false);
        }, 2000);
      })
      .catch((err) => {
        console.error("Failed to copy primers:", err);
      });
  }, [result]);

  function formatPrimers(primer) {
    let seq = formatPrimerSequence(primer);
    return <span className="primer-sequence">{seq}</span>;
  }

  return (
    <div className="sequence-results">
      {/* Messages */}
      {result.messages && result.messages.length > 0 && (
        <div className="messages">
          {result.messages.map((msg, index) => (
            <div key={index}>{msg}</div>
          ))}
        </div>
      )}
      {/* Display MTK Part Numbers */}
      <div className="mtk-part-info">
        {mtkPartLeft === mtkPartRight ? (
          <p>
            <strong>MTK Part Number:</strong> {mtkPartLeft}
          </p>
        ) : (
          <p>
            <strong>MTK Part Number Left:</strong> {mtkPartLeft} <br />
            <strong>MTK Part Number Right:</strong> {mtkPartRight}
          </p>
        )}
      </div>
      {/* Restriction Sites */}
      {result.restriction_sites && result.restriction_sites.length > 0 && (
        <RestrictionSiteSummary sites={result.restriction_sites} />
      )}

      {/* PCR Reactions */}
      {result.PCR_reactions && Object.keys(result.PCR_reactions).length > 0 && (
        <div className="pcr-summary">
          <div className="section-header">
            <h4>PCR Reactions</h4>
            {/* New Button to Copy Primers */}
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
                  ([reaction, primers], index) => (
                    <tr key={index}>
                      <td>{reaction}</td>
                      <td className="primer-cell">
                        {formatPrimers(primers.forward)}
                      </td>
                      <td className="primer-cell">
                        {formatPrimers(primers.reverse)}
                      </td>
                    </tr>
                  )
                )}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Errors */}
      {result.errors && (
        <div className="error-message">
          <strong>Error:</strong> {result.errors}
        </div>
      )}
    </div>
  );
}

export default ResultTab;
