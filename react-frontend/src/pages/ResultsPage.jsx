// src/pages/ResultsPage.jsx
import React, { useState, useEffect, useCallback } from "react";
import { useNavigate } from "react-router-dom";
import Results from "../components/Results/Results";
import { monitorProtocolProgress } from "../api/api";

function ResultsPage({ results }) {
  const [localResults, setLocalResults] = useState(results);
  // progress is an object keyed by sequence id (or a global key)
  const [progress, setProgress] = useState({});
  const navigate = useNavigate();

  useEffect(() => {
    // If no results passed via props, try sessionStorage.
    if (!results) {
      const storedResults = sessionStorage.getItem("results");
      if (storedResults) {
        setLocalResults(JSON.parse(storedResults));
      } else {
        // If we truly have nothing, redirect to the form.
        navigate("/");
      }
    }
  }, [results, navigate]);

  // Use useCallback to ensure a stable reference for the callback.
  const onStatusUpdate = useCallback((statusData) => {
    // Update per-sequence progress if sequenceId is present, otherwise update global progress.
    if (statusData.sequenceId !== undefined) {
      setProgress((prev) => ({
        ...prev,
        [statusData.sequenceId]: statusData,
      }));
    } else {
      setProgress((prev) => ({
        ...prev,
        global: statusData,
      }));
    }
  }, []);

  // Set up SSE connection if a jobId is stored (assuming FormPage stored it).
  useEffect(() => {
    const jobId = sessionStorage.getItem("jobId");
    if (jobId) {
      const eventSource = monitorProtocolProgress(jobId, onStatusUpdate);
      // Clean up the event source on unmount.
      return () => {
        eventSource.close();
      };
    }
  }, [onStatusUpdate]);

  return (
    <div className="output-container">
      {/* Display a global progress indicator if available */}
      {progress.global && (
        <div className="global-progress">
          <p>
            Global Progress: {progress.global.percentage}% —{" "}
            {progress.global.message}
          </p>
        </div>
      )}
      {localResults ? (
        <Results data={localResults} progress={progress} />
      ) : (
        <p className="initialization-message">Loading...</p>
      )}
    </div>
  );
}

export default ResultsPage;
