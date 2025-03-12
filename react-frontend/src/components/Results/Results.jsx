// import React, { useState } from "react";
// import PrimerAnatomy from "./PrimerAnatomy";
import ResultTabs from "./ResultTabs";
import "../../styles/Results.css";

function Results({ data }) {
  // const [showPrimerAnatomy, setShowPrimerAnatomy] = useState(false);

  if (!data || Object.keys(data).length === 0) return null;

  // Convert data object to array for use with ResultTabs
  const resultsArray = Object.entries(data).map(
    ([sequenceNumber, sequenceData]) => ({
      sequenceNumber,
      ...sequenceData,
    })
  );

  return (
    <div className="protocol-results">
      <h2>Golden Gate Protocol Results</h2>

      {/* Toggle Primer Anatomy View
      <div className="result-actions">
        <button
          className="btn btn-secondary"
          onClick={() => setShowPrimerAnatomy(!showPrimerAnatomy)}
        >
          {showPrimerAnatomy ? "Hide" : "Show"} Primer Anatomy
        </button>
      </div> */}

      <ResultTabs results={resultsArray} />

      {/* {showPrimerAnatomy && (
        <div className="primer-anatomy-container">
          <PrimerAnatomy />
        </div>
      )} */}
    </div>
  );
}

export default Results;
