// src/components/Results/Results.jsx
import ResultTabs from "./ResultTabs";
import "../../styles/Results.css";

function Results({ data, progress }) {

  if (!data || Object.keys(data).length === 0) return null;

  // Convert the data object to an array and attach progress info (if any)
  const resultsArray = Object.entries(data).map(
    ([sequenceKey, sequenceData]) => ({
      sequenceKey,
      ...sequenceData,
      progress: progress ? progress[sequenceKey] : null,
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
