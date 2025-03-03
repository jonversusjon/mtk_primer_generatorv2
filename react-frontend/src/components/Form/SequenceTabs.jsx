import React, { useState } from 'react';
import SequenceTab from './SequenceTab';
import '../../styles/SequenceTabs.css';

const MTK_PART_NUMS = ['', '1', '2', '3', '3a', '3b', '3c', '3d', '3e', '3f',
                       '3g', '4', '4a', '4b', '4aII', '3bI', '3bII', '5', '6', '7', '8', '8a', '8b'];

function SequenceTabs({ sequences, updateSequence, addSequence, removeSequence }) {
  const [activeTab, setActiveTab] = useState(0);
  
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
            disabled={sequences.length <= 1}
          >
            –
          </button>
          <button 
            type="button" 
            className="btn btn-sm btn-outline-success"
            onClick={addSequence}
            disabled={sequences.length >= 10}
          >
            +
          </button>
        </div>
      </div>
      
      {/* Sequence Tabs Container */}
      <div className="sequence-tabs-container">
        {/* Navigation Buttons */}
        <div className="tab-buttons">
          {sequences.map((_, index) => (
            <button
              key={index}
              type="button"
              className={`tab-button ${activeTab === index ? 'active' : ''}`}
              onClick={() => setActiveTab(index)}
            >
              {index + 1}
            </button>
          ))}
        </div>
        
        {/* Tab Contents */}
        <div className="tab-content">
          {sequences.map((sequence, index) => (
            <div
              key={index}
              className={`tab-pane ${activeTab === index ? 'active' : ''}`}
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