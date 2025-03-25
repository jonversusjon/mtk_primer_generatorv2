import React, { useCallback, useState } from "react";

// Helper to format primer sequences consistently.
const formatPrimerSequence = (primer) => {
  if (!primer) return "None";
  if (typeof primer === "string") return primer;
  if (Array.isArray(primer)) return primer.join(", ");
  if (primer.sequence) return primer.sequence;
  return "None";
};


const PcrReaction = ({ pcrReactions }) => {
  const [copied, setCopied] = useState(false);

  const copyPrimersToClipboard = useCallback(() => {
    if (!pcrReactions) return;

    const rows = Object.entries(pcrReactions).flatMap(([reactionName, primers]) => {
      const forwardSeq = formatPrimerSequence(primers.forward);
      const reverseSeq = formatPrimerSequence(primers.reverse);
      return [
        `${reactionName}_FWD\t${forwardSeq}`,
        `${reactionName}_REV\t${reverseSeq}`,
      ];
    });
    const finalText = rows.join("\n");

    navigator.clipboard
      .writeText(finalText)
      .then(() => {
        setCopied(true);
        setTimeout(() => setCopied(false), 2000);
      })
      .catch((err) => console.error("Failed to copy primers:", err));
  }, [pcrReactions]);

  if (!pcrReactions || Object.keys(pcrReactions).length === 0) {
    return null;
  }

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
            {Object.entries(pcrReactions).map(([reaction, primers], idx) => (
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
};

export default PcrReaction;
