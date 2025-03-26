// ResultTabs.jsx
import React, { useState, useEffect } from "react";
import ResultTab from "./ResultTab";

const ResultTabs = () => {
  // Load initial sequences from sessionStorage
  const [results, setResults] = useState([]);

  useEffect(() => {
    const savedFormData = sessionStorage.getItem("formData");
    if (savedFormData) {
      try {
        const parsed = JSON.parse(savedFormData);
        if (parsed.sequencesToDomesticate?.length > 0) {
          const initialResults = parsed.sequencesToDomesticate.map(
            (seq, i) => ({
              id: i,
              sequence: seq.sequence,
              primerName: seq.primerName || `Sequence ${i + 1}`,
              placeholder: true,
            })
          );
          setResults(initialResults);
        }
      } catch (error) {
        console.error("Error parsing formData from sessionStorage:", error);
      }
    }
  }, []);

  const [activeTab, setActiveTab] = useState(0);

  // If activeTab is out of bounds, adjust it.
  useEffect(() => {
    if (activeTab >= results.length && results.length > 0) {
      setActiveTab(results.length - 1);
    }
  }, [results.length, activeTab]);

  // If no results available yet, show a loading message.
  if (results.length === 0) {
    return <p className="initialization-message">Loading...</p>;
  }

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
              className={`tab-button results-tab-button ${
                isActive ? "active" : ""
              }`}
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
            {/* Each ResultTab receives its own sequence index so it subscribes to the proper SSE channel */}
            <ResultTab result={result} sequenceIdx={index} />
          </div>
        ))}
      </div>
    </div>
  );
};

export default ResultTabs;
