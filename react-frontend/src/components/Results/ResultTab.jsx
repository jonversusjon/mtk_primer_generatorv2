import React from "react";
import RestrictionSiteSummary from "./RestrictionSiteSummary";

// Helper to format primer sequences with a special span for styling
function formatPrimers(primer) {
  let seq = "";
  if (!primer) seq = "None";
  else if (typeof primer === "string") seq = primer;
  else if (Array.isArray(primer)) seq = primer.join(", ");
  else if (primer.sequence) seq = primer.sequence;
  else seq = "None";
  return <span className="primer-sequence">{seq}</span>;
}

function ResultTab({ result }) {
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

      {/* Restriction Sites */}
      {result.restriction_sites && result.restriction_sites.length > 0 && (
        <RestrictionSiteSummary sites={result.restriction_sites} />
      )}

      {/* Mutations */}
      {result.mutations &&
        result.mutations.all_mutation_options &&
        result.mutations.all_mutation_options.length > 0 && (
          <div className="mutations-summary">
            <h4>Mutations</h4>
            <div className="table-container">
              <table>
                <thead>
                  <tr>
                    <th>Position</th>
                    <th>Original</th>
                    <th>Mutated</th>
                    <th>Type</th>
                  </tr>
                </thead>
                <tbody>
                  {result.mutations.all_mutation_options.map(
                    (mutation, index) => (
                      <tr key={index}>
                        <td>{mutation.position}</td>
                        <td>{mutation.original_sequence}</td>
                        <td>{mutation.mutated_sequence}</td>
                        <td>{mutation.type}</td>
                      </tr>
                    )
                  )}
                </tbody>
              </table>
            </div>
          </div>
        )}

      {/* PCR Reactions */}
      {result.PCR_reactions && Object.keys(result.PCR_reactions).length > 0 && (
        <div className="pcr-summary">
          <h4>PCR Reactions</h4>
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
