// ResultTabs.jsx
import React, { useState, useEffect } from "react";
import ResultTab from "./ResultTab";

function ResultTabs({ results }) {
  const [activeTab, setActiveTab] = useState(0);

  // Ensure active tab is valid after changes in results
  useEffect(() => {
    if (activeTab >= results.length && results.length > 0) {
      setActiveTab(results.length - 1);
    }
  }, [results.length, activeTab]);

  return (
    <div className="results-section">
      <div className="section-header">
        <div className="section-title">
        </div>
      </div>

      <div className="result-tabs-container">
        {/* Tab Navigation */}
        <div className="tab-buttons">
          {results.map((result, index) => (
            <button
              key={index}
              type="button"
              className={`tab-button results-tab-button ${
                activeTab === index ? "active" : ""
              }`}
              onClick={() => setActiveTab(index)}
              role="tab"
              aria-selected={activeTab === index}
              aria-controls={`result-tab-${index}`}
              id={`result-tab-button-${index}`}
            >
              {result.primerName?.trim()
                ? result.primerName
                : `Sequence ${index + 1}`}
              {/* Show a progress badge if processing is still underway */}
              {result.progress && result.progress.percentage < 100 && (
                <span className="progress-badge">
                  {result.progress.percentage}%
                </span>
              )}
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

export default ResultTabs;
