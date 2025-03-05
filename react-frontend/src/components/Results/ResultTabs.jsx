import React, { useState, useEffect } from "react";

function ResultTabs({ results, removeResult, addResult }) {
  const [activeTab, setActiveTab] = useState(0);

  // Ensure active tab is valid after removing results
  useEffect(() => {
    if (activeTab >= results.length && results.length > 0) {
      setActiveTab(results.length - 1);
    }
  }, [results.length, activeTab]);

  return (
    <div className="results-section">
      <div className="section-header">
        <div className="section-title">
          <h3>Results</h3>
        </div>
      </div>

      <div className="result-tabs-container">
        {/* Tab Navigation */}
        <div className="tab-buttons">
          {results.map((_, index) => (
            <button
              key={index}
              type="button"
              className={`tab-button ${activeTab === index ? "active" : ""}`}
              onClick={() => setActiveTab(index)}
              role="tab"
              aria-selected={activeTab === index}
              aria-controls={`result-tab-${index}`}
              id={`result-tab-button-${index}`}
            >
              {index + 1}
            </button>
          ))}
        </div>

        {/* Tab Content Panels */}
        <div className="tab-content">
          {results.map((result, index) => (
            <div
              key={index}
              className={`tab-pane ${activeTab === index ? "active" : ""}`}
              id={`result-tab-${index}`}
              role="tabpanel"
              aria-labelledby={`result-tab-button-${index}`}
              hidden={activeTab !== index}
            >
              <ResultTab result={result} index={index} />
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}

function ResultTab({ result, index }) {
  return (
    <div className="result-content">
      <h4>Result {index + 1}</h4>
      <pre>{JSON.stringify(result, null, 2)}</pre>
    </div>
  );
}

export default ResultTabs;
