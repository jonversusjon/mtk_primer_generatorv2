// Results/ResultTabs.jsx
import React, { useState, useEffect } from "react";
import ResultTab from "./ResultTab"; // Ensure this path is correct

const ResultTabs = ({ jobId }) => {
  // Load sequences from sessionStorage (each sequence gets an id that serves as sequenceIdx)
  const [sequences, setSequences] = useState([]);
  const [activeTab, setActiveTab] = useState(0); // Keep track of the active tab index

  useEffect(() => {
    // Attempt to load sequence information from sessionStorage
    const savedFormData = sessionStorage.getItem("formData");
    if (savedFormData) {
      try {
        const parsed = JSON.parse(savedFormData);
        // Assuming sequencesToDomesticate holds the array of sequence objects
        if (parsed.sequencesToDomesticate?.length > 0) {
          // Map sequences to include an 'id' (which corresponds to sequenceIdx) and a display name
          const seqs = parsed.sequencesToDomesticate.map((seq, i) => ({
            id: i, // Use the index as the unique ID and sequenceIdx
            // Use primerName if available, otherwise fallback to a generic name
            displayName: seq.primerName?.trim() || `Sequence ${i + 1}`,
          }));
          setSequences(seqs);
        } else {
          console.warn("No sequences found in sessionStorage formData.");
          setSequences([]); // Ensure sequences is an empty array if none are found
        }
      } catch (error) {
        console.error("Error parsing formData from sessionStorage:", error);
        setSequences([]); // Reset sequences on error
      }
    } else {
      console.warn("No formData found in sessionStorage.");
      setSequences([]); // Ensure sequences is an empty array if no data found
    }
  }, []); // Run this effect only once on component mount

  // Adjust activeTab if it becomes out of bounds (e.g., if sequences load async or change)
  useEffect(() => {
    if (sequences.length > 0) {
      // If activeTab is invalid, reset to 0 or the last valid index
      if (activeTab >= sequences.length || activeTab < 0) {
        setActiveTab(0); // Default to the first tab if current is invalid
      }
    } else {
      // If there are no sequences, reset activeTab to 0
      setActiveTab(0);
    }
  }, [sequences, activeTab]);

  // Add debug logging for state changes
  useEffect(() => {
    console.log("Active tab changed to:", activeTab);
  }, [activeTab]);

  // Display a loading message or placeholder if sequences haven't loaded yet
  if (sequences.length === 0) {
    // Added check for jobId as well, perhaps loading is dependent on it too
    return <p className="p-4 text-gray-500 dark:text-gray-400">Loading sequence data{jobId ? ` for Job ${jobId}` : ''}...</p>;
  }

  return (
    <div className="results-section w-full max-w-4xl mx-auto"> 
      {/* Tab Buttons Container */}
      <div className="flex border-b border-gray-300 dark:border-gray-700 mb-4 overflow-x-auto">
        {sequences.map((seq, index) => {
          const isActive = activeTab === index;
          return (
            <button
              key={seq.id}
              type="button"
              role="tab"
              className={`py-2 px-4 text-sm font-medium text-center border-b-2 whitespace-nowrap ${
                isActive
                  ? "border-blue-500 text-blue-600 dark:text-blue-400 dark:border-blue-400"
                  : "border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300 dark:text-gray-400 dark:hover:text-gray-300 dark:hover:border-gray-600"
              }`}
              onClick={() => {
                console.log("Setting activeTab to:", index);
                setActiveTab(index);
              }}
              aria-selected={isActive}
              aria-controls={`tab-content-${seq.id}`}
              id={`tab-button-${seq.id}`}
            >
              {seq.displayName}
            </button>
          );
        })}
      </div>

      {/* Tab Content Area */}
      <div className="tab-content relative">
        {/* Render ALL ResultTab components, but control visibility with inline styles instead of hidden attribute */}
        {sequences.map((seq, index) => (
          <div
            key={seq.id}
            role="tabpanel"
            aria-labelledby={`tab-button-${seq.id}`}
            id={`tab-content-${seq.id}`}
            className="tab-pane"
            style={{ 
              display: activeTab === index ? 'block' : 'none',
              // Optional debugging border to verify visibility
              border: activeTab === index ? '2px solid blue' : 'none'
            }}
          >
            <ResultTab jobId={jobId} sequenceIdx={seq.id} />
          </div>
        ))}
      </div>
    </div>
  );
};

export default ResultTabs;