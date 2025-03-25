import React, { useState, useEffect } from "react";
import ResultTab from "./ResultTab";

const ResultTabs = ({ results }) => {
  const [activeTab, setActiveTab] = useState(0);

  useEffect(() => {
    if (activeTab >= results.length && results.length > 0) {
      setActiveTab(results.length - 1);
    }
  }, [results.length, activeTab]);

  return (
    <div className="results-section">
      <div className="tab-buttons">
        {results.map((result, index) => {
          const isActive = activeTab === index;
          const tabLabel = result.primerName?.trim() || `Sequence ${index + 1}`;
          return (
            <button
            key={index}
            type="button"
            role="tab"
            className={`tab-button results-tab-button ${isActive ? "active" : ""}`}
            onClick={() => setActiveTab(index)}
            aria-selected={isActive}
            aria-controls={`tab-content-${index}`}
            id={`tab-button-${index}`}
          >
              {tabLabel}
            </button>
          );
        })}
      </div>
      <div className="tab-content">
        {results.map((result, index) => (
          <div
            key={index}
            className={`tab-pane ${activeTab === index ? "active" : ""}`}
            role="tabpanel"
            hidden={activeTab !== index}
          >
            {/* Pass jobId and sequence index to ResultTab */}
            <ResultTab result={result} sequenceIdx={index} />
          </div>
        ))}
      </div>
    </div>
  );
};

export default ResultTabs;
