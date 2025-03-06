import React, { useState, useEffect } from "react";
import RestrictionSiteSummary from "./RestrictionSiteSummary";
import PrimerAnatomy from "./PrimerAnatomy";
import "../../styles/Results.css";

function Results({ data }) {
  const [activeTab, setActiveTab] = useState(0);
  const [showPrimerAnatomy, setShowPrimerAnatomy] = useState(false);

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

      {/* Toggle Primer Anatomy View */}
      <div className="result-actions">
        <button
          className="btn btn-secondary"
          onClick={() => setShowPrimerAnatomy(!showPrimerAnatomy)}
        >
          {showPrimerAnatomy ? "Hide" : "Show"} Primer Anatomy
        </button>
      </div>

      {showPrimerAnatomy && (
        <div className="primer-anatomy-container">
          <PrimerAnatomy />
        </div>
      )}

      <ResultTabs
        results={resultsArray}
        activeTab={activeTab}
        setActiveTab={setActiveTab}
        removeResult={() => {}}
        addResult={() => {}}
      />
    </div>
  );
}

function ResultTabs({
  results,
  activeTab,
  setActiveTab,
  removeResult,
  addResult,
}) {
  // Ensure active tab is valid after removing results
  useEffect(() => {
    if (activeTab >= results.length && results.length > 0) {
      setActiveTab(results.length - 1);
    }
  }, [results.length, activeTab, setActiveTab]);

  return (
    <div className="results-section">

      <div className="result-tabs-container">
        {/* Tab Navigation */}
        <div className="tab-buttons">
          {results.map((result, index) => (
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
              Sequence {result.sequenceNumber}
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
              <ResultTab result={result} />
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}

function ResultTab({ result }) {
  return (
    <div className="sequence-results">
      {/* Messages */}
      {result.messages && result.messages.length > 0 && (
        <div className="messages">
          {result.messages.map((msg, index) => (
            <div key={index}>{msg}</div>
          ))}
        </div>
      )}

      {/* Restriction Sites */}
      {result.restriction_sites && result.restriction_sites.length > 0 && (
        <RestrictionSiteSummary sites={result.restriction_sites} />
      )}

      {/* Mutations */}
      {result.mutations &&
        result.mutations.all_mutation_options &&
        result.mutations.all_mutation_options.length > 0 && (
          <div className="mutations-summary">
            <h4>Mutations</h4>
            <table>
              <thead>
                <tr>
                  <th>Position</th>
                  <th>Original</th>
                  <th>Mutated</th>
                  <th>Type</th>
                </tr>
              </thead>
              <tbody>
                {result.mutations.all_mutation_options.map(
                  (mutation, index) => (
                    <tr key={index}>
                      <td>{mutation.position}</td>
                      <td>{mutation.original_sequence}</td>
                      <td>{mutation.mutated_sequence}</td>
                      <td>{mutation.type}</td>
                    </tr>
                  )
                )}
              </tbody>
            </table>
          </div>
        )}

      {/* PCR Reactions */}
      {result.PCR_reactions && Object.keys(result.PCR_reactions).length > 0 && (
        <div className="pcr-summary">
          <h4>PCR Reactions</h4>
          <table>
            <thead>
              <tr>
                <th>Reaction</th>
                <th>Mutation Primers</th>
                <th>Edge Primers</th>
              </tr>
            </thead>
            <tbody>
              {Object.entries(result.PCR_reactions).map(
                ([reaction, primers], index) => (
                  <tr key={index}>
                    <td>{reaction}</td>
                    <td>
                      {primers.mutation_primers
                        ? primers.mutation_primers.join(", ")
                        : "None"}
                    </td>
                    <td>
                      {primers.edge_primers
                        ? primers.edge_primers.join(", ")
                        : "None"}
                    </td>
                  </tr>
                )
              )}
            </tbody>
          </table>
        </div>
      )}

      {/* Errors */}
      {result.errors && (
        <div className="error-message">
          <strong>Error:</strong> {result.errors}
        </div>
      )}
    </div>
  );
}

export default Results;
