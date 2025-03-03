import React, { useState } from 'react';
import SequenceTab from './SequenceTab';
import '../../styles/SequenceTabs.css';

const MTK_PART_NUMS = ['', '1', '2', '3', '3a', '3b', '3c', '3d', '3e', '3f',
                       '3g', '4', '4a', '4b', '4aII', '3bI', '3bII', '5', '6', '7', '8', '8a', '8b'];

function SequenceTabs({ sequences, updateSequence, addSequence, removeSequence }) {
  const [activeTab, setActiveTab] = useState(0);
  
  return (
    <div className="form-sub-container">
      <div className="section-header">
        <h2 className="form-sub-container-title">Sequences to Domesticate</h2>
        <div id="sequence-control-buttons">
          <button
            type="button"
            className="sequence-control-btn"
            onClick={removeSequence}
            disabled={sequences.length <= 1}
          >
            –
          </button>
          <button
            type="button"
            className="sequence-control-btn"
            onClick={addSequence}
            disabled={sequences.length >= 10}
          >
            +
          </button>
        </div>
      </div>
      
      {/* Sequence Tabs Container */}
      <div id="sequence-inputs-container" className="tabs-container">
        {/* Navigation Buttons */}
        <div className="sequence-tabs-nav">
          {sequences.map((_, index) => (
            <button
              key={index}
              type="button"
              className={`sequence-tab-btn nav-narrow ${activeTab === index ? 'active' : ''}`}
              onClick={() => setActiveTab(index)}
            >
              {index + 1}
            </button>
          ))}
        </div>
        
        {/* Tab Contents */}
        <div className="sequence-tabs-content">
          {sequences.map((sequence, index) => (
            <div
              key={index}
              className={`sequence-tab-content ${activeTab === index ? 'active' : 'hidden'}`}
            >
              <SequenceTab
                sequence={sequence}
                index={index}
                updateSequence={(field, value) => updateSequence(index, field, value)}
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