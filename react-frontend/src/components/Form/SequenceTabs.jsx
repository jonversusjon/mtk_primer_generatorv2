import React, { useState } from "react";
import SequenceTab from "./SequenceTab";
import "../../styles/SequenceTabs.css";

// Constants should be outside the component to avoid recreating on each render
const MTK_PART_NUMS = [
  "",
  "1",
  "2",
  "3",
  "3a",
  "3b",
  "3c",
  "3d",
  "3e",
  "3f",
  "3g",
  "4",
  "4a",
  "4b",
  "4aII",
  "3bI",
  "3bII",
  "5",
  "6",
  "7",
  "8",
  "8a",
  "8b",
];

function SequenceTabs({
  sequencesToDomesticate,
  updateSequence,
  addSequence,
  removeSequence,
}) {
  const [activeTab, setActiveTab] = useState(0);

  // Ensure active tab is valid after removing sequences
  React.useEffect(() => {
    if (
      activeTab >= sequencesToDomesticate.length &&
      sequencesToDomesticate.length > 0
    ) {
      setActiveTab(sequencesToDomesticate.length - 1);
    }
  }, [sequencesToDomesticate.length, activeTab]);

  return (
    <div className="sequences-section">
      <div className="section-header">
        <div className="section-title">
          <h3>Sequences to Domesticate</h3>
        </div>

        <div className="sequence-controls">
          <button
            type="button"
            className="btn btn-sm btn-outline-danger"
            onClick={removeSequence}
            disabled={sequencesToDomesticate.length <= 1}
            aria-label="Remove sequence"
          >
            –
          </button>
          <button
            type="button"
            className="btn btn-sm btn-outline-success"
            onClick={addSequence}
            disabled={sequencesToDomesticate.length >= 10}
            aria-label="Add sequence"
          >
            +
          </button>
        </div>
      </div>

      <div className="sequence-tabs-container">
        {/* Tab Navigation */}
        <div className="tab-buttons">
          {sequencesToDomesticate.map((_, index) => (
            <button
              key={index}
              type="button"
              className={`tab-button ${activeTab === index ? "active" : ""}`}
              onClick={() => setActiveTab(index)}
              role="tab"
              aria-selected={activeTab === index}
              aria-controls={`sequence-tab-${index}`}
              id={`sequence-tab-button-${index}`}
            >
              {index + 1}
            </button>
          ))}
        </div>

        {/* Tab Content Panels */}
        <div className="tab-content">
          {sequencesToDomesticate.map((sequence, index) => (
            <div
              key={index}
              className={`tab-pane ${activeTab === index ? "active" : ""}`}
              id={`sequence-tab-${index}`}
              role="tabpanel"
              aria-labelledby={`sequence-tab-button-${index}`}
              hidden={activeTab !== index}
            >
              <SequenceTab
                sequence={sequence}
                index={index}
                updateSequence={updateSequence(index)}
                mtkPartOptions={MTK_PART_NUMS}
              />
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}

export default SequenceTabs;
