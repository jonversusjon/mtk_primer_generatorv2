import React from "react";
import ResultTabs from "./ResultTabs";
import "../../styles/Results.css";

const Results = ({ data }) => {
  if (!data || (Array.isArray(data) && data.length === 0)) return null;

  // If data is an object, convert it to an array.
  const resultsArray = Array.isArray(data) ? data : Object.values(data);

  return (
    <div className="protocol-results">
      <h2>Golden Gate Protocol Results</h2>
      <ResultTabs results={resultsArray} />
    </div>
  );
};

export default Results;
