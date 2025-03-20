import React, { useState, useEffect } from "react";
import ResultTab from "./ResultTab";

/*
  ResultTabs renders the tab navigation and content panels.
  - The tab label is updated to show "(Processing)" if the result is a placeholder
    or its progress percentage is less than 100.
  - It also ensures the active tab index remains valid when the results change.
*/
const ResultTabs = ({ results }) => {
  const [activeTab, setActiveTab] = useState(0);

  // Update active tab index if the number of results changes.
  useEffect(() => {
    if (activeTab >= results.length && results.length > 0) {
      setActiveTab(results.length - 1);
    }
  }, [results.length, activeTab]);

  return (
    <div className="results-section">
      <div className="section-header">
        <div className="section-title"></div>
      </div>
      <div className="result-tabs-container">
        {/* Tab Navigation */}
        <div className="tab-buttons">
          {results.map((result, index) => {
            // Build tab label:
            let tabLabel = result.primerName?.trim()
              ? result.primerName
              : `Sequence ${index + 1}`;
            // Append a processing indicator if placeholder or not yet complete.
            if (result.placeholder || (result.progress && result.progress.percentage < 100)) {
              tabLabel += " (Processing)";
            }
            return (
              <button
                key={index}
                type="button"
                className={`tab-button results-tab-button ${activeTab === index ? "active" : ""}`}
                onClick={() => setActiveTab(index)}
                role="tab"
                aria-selected={activeTab === index}
                aria-controls={`result-tab-${index}`}
                id={`result-tab-button-${index}`}
              >
                {tabLabel}
                {result.progress && result.progress.percentage < 100 && (
                  <span className="progress-badge">{result.progress.percentage}%</span>
                )}
              </button>
            );
          })}
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
              <ResultTab result={result} />
            </div>
          ))}
        </div>
      </div>
    </div>
  );
};

export default ResultTabs;
