import React, { useState } from "react";
import RestrictionSiteSummary from "./RestrictionSiteSummary";
import PrimerAnatomy from "./PrimerAnatomy";
import "../../styles/Results.css";

function Results({ data }) {
  const [showPrimerAnatomy, setShowPrimerAnatomy] = useState(false);

  if (!data) return null;

  const { restriction_sites, primers, mutations, sequence_analysis } = data;

  return (
    <div className="protocol-results">
      <h2>Golden Gate Protocol Results</h2>

      {/* Toggle Primer Anatomy View */}
      <div className="result-actions">
        <button
          className="btn btn-secondary"
          onClick={() => setShowPrimerAnatomy(!showPrimerAnatomy)}
        >
          {showPrimerAnatomy ? "Hide" : "Show"} Primer Anatomy
        </button>
      </div>

      {/* Show Primer Anatomy if toggled */}
      {showPrimerAnatomy && (
        <div className="primer-anatomy-container">
          <PrimerAnatomy />
        </div>
      )}

      {/* Show Restriction Site Summary if available */}
      {restriction_sites && restriction_sites.length > 0 && (
        <RestrictionSiteSummary sites={restriction_sites} />
      )}

      {/* Show Primers if available */}
      {primers && primers.length > 0 && (
        <div className="primers-summary">
          <h3>Designed Primers</h3>
          <table>
            <thead>
              <tr>
                <th>Name</th>
                <th>Sequence</th>
                <th>Length</th>
              </tr>
            </thead>
            <tbody>
              {primers.map((primer, index) => (
                <tr key={index}>
                  <td>{primer[0]}</td>
                  <td className="primer-sequence">{primer[1]}</td>
                  <td>{primer[1].length} bp</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      {/* Show Mutations if available */}
      {mutations && mutations.length > 0 && (
        <div className="mutations-summary">
          <h3>Mutations</h3>
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
              {mutations.map((mutation, index) => (
                <tr key={index}>
                  <td>{mutation.position}</td>
                  <td>{mutation.original}</td>
                  <td>{mutation.mutated}</td>
                  <td>{mutation.type}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      {/* Show any additional analysis data */}
      {sequence_analysis && sequence_analysis.length > 0 && (
        <div className="sequence-analysis-summary">
          <h3>Sequence Analysis</h3>
          <div className="analysis-details">
            {sequence_analysis.map((analysis, index) => (
              <div key={index} className="analysis-item">
                <h4>Sequence {index + 1}</h4>
                <p>
                  <strong>Processed Length:</strong>{" "}
                  {analysis.processed_sequence?.length || 0} bp
                </p>
                {analysis.restriction_sites &&
                  analysis.restriction_sites.length > 0 && (
                    <p>
                      <strong>Restriction Sites:</strong>{" "}
                      {analysis.restriction_sites.length}
                    </p>
                  )}
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}

export default Results;
