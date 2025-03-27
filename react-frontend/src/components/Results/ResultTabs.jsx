// ResultTabs.jsx
import React, { useState, useEffect } from "react";
import ResultTab from "./ResultTab";

const ResultTabs = ({ jobId }) => {
  // Load sequences to domesticate from sessionStorage
  const [sequences, setSequences] = useState([]);

  useEffect(() => {
    const savedFormData = sessionStorage.getItem("formData");
    if (savedFormData) {
      try {
        const parsed = JSON.parse(savedFormData);
        if (parsed.sequencesToDomesticate?.length > 0) {
          // Map each sequence to an object with an id (which serves as sequenceIdx)
          const seqs = parsed.sequencesToDomesticate.map((seq, i) => ({
            id: i,
            primerName: seq.primerName || `Sequence ${i + 1}`,
          }));
          setSequences(seqs);
        }
      } catch (error) {
        console.error("Error parsing formData from sessionStorage:", error);
      }
    }
  }, []);

  const [activeTab, setActiveTab] = useState(0);

  // Ensure activeTab is within bounds when sequences update.
  useEffect(() => {
    if (activeTab >= sequences.length && sequences.length > 0) {
      setActiveTab(sequences.length - 1);
    }
  }, [sequences, activeTab]);

  if (sequences.length === 0) {
    return <p className="initialization-message">Loading...</p>;
  }

  return (
    <div className="results-section">
      <div className="tab-buttons">
        {sequences.map((seq, index) => {
          const isActive = activeTab === index;
          const tabLabel = seq.primerName?.trim() || `Sequence ${index + 1}`;
          return (
            <button
              key={seq.id}
              type="button"
              role="tab"
              className={`tab-button results-tab-button ${isActive ? "active" : ""}`}
              onClick={() => setActiveTab(index)}
              aria-selected={isActive}
              aria-controls={`tab-content-${seq.id}`}
              id={`tab-button-${seq.id}`}
            >
              {tabLabel}
            </button>
          );
        })}
      </div>
      <div className="tab-content">
        {sequences.map((seq, index) => (
          <div
            key={seq.id}
            className={`tab-pane ${activeTab === index ? "active" : ""}`}
            role="tabpanel"
            hidden={activeTab !== index}
          >
            {/* Pass jobId and sequenceIdx to each ResultTab */}
            <ResultTab jobId={jobId} sequenceIdx={seq.id} />
          </div>
        ))}
      </div>
    </div>
  );
};

export default ResultTabs;
