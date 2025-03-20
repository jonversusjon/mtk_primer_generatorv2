import React from "react";
import ResultTabs from "./ResultTabs";
import "../../styles/Results.css";

/*
  This component receives `data` and `progress` props.
  - If `data` is an array (placeholders from sessionStorage) or an object (final data), we map it
    into a consistent results array.
  - For each sequence, we attach any progress info from SSE using the sequence's key.
*/
const Results = ({ data, progress }) => {
  // Check for empty data (array or object)
  if (
    !data ||
    (Array.isArray(data) && data.length === 0) ||
    (!Array.isArray(data) && Object.keys(data).length === 0)
  )
    return null;

  let resultsArray;
  if (Array.isArray(data)) {
    // Data is already an array (placeholders)
    resultsArray = data.map((item, index) => ({
      sequenceKey: item.id || index,
      ...item,
      progress: progress ? progress[item.id || index] : null,
    }));
  } else {
    // Data is an object (final data keyed by sequence)
    resultsArray = Object.entries(data).map(([sequenceKey, sequenceData]) => ({
      sequenceKey,
      ...sequenceData,
      progress: progress ? progress[sequenceKey] : null,
    }));
  }

  return (
    <div className="protocol-results">
      <h2>Golden Gate Protocol Results</h2>
      <ResultTabs results={resultsArray} />
    </div>
  );
};

export default Results;
