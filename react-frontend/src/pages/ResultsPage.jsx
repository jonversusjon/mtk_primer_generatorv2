// ResultsPage.jsx
import React, { useState, useEffect, useCallback, useMemo } from "react";
import { useNavigate } from "react-router-dom";
import Results from "../components/Results/Results";
import { monitorProtocolProgress } from "../api/api";

function ResultsPage({ results }) {
  const navigate = useNavigate();

  // Load final results (if any) from sessionStorage.
  const [finalResults, setFinalResults] = useState(() => {
    const savedResults = sessionStorage.getItem("results");
    return savedResults ? JSON.parse(savedResults) : null;
  });

  // Load the initial message from sessionStorage (e.g., "Primer design started")
  const [initialMessage] = useState(() => {
    return sessionStorage.getItem("initialMessage") || "Primer design started...";
  });

  // Build placeholder objects from formData in sessionStorage.
  const placeholders = useMemo(() => {
    const savedFormData = sessionStorage.getItem("formData");
    if (savedFormData) {
      const parsedFormData = JSON.parse(savedFormData);
      if (parsedFormData.sequencesToDomesticate) {
        return parsedFormData.sequencesToDomesticate.map((seq, index) => ({
          id: index,
          placeholder: true,
          sequence: seq.sequence,
          primerName: seq.primerName || `Sequence ${index + 1}`,
        }));
      }
    }
    return [];
  }, []);

  // Local progress state for SSE updates.
  // This will store global progress as well as per-sequence progress (if available).
  const [progress, setProgress] = useState({});

  // Determine what data to display:
  // If final results are available, they take precedence.
  // Otherwise, use the results passed as a prop, or fall back to the placeholders.
  const dataToDisplay = finalResults || results || placeholders;

  // Redirect to the form if no data is available.
  useEffect(() => {
    if (!dataToDisplay || dataToDisplay.length === 0) {
      navigate("/");
    }
  }, [dataToDisplay, navigate]);

  // Callback to handle SSE updates.
  // Updates the progress state and, if a global update reaches 100%,
  // updates the final results.
  const onStatusUpdate = useCallback((statusData) => {
    // If the status data contains a sequenceId, use it; otherwise, it's global.
    const key = statusData.sequenceId !== undefined ? statusData.sequenceId : "global";
    setProgress((prev) => ({ ...prev, [key]: statusData }));

    // Update final results if global progress reaches 100%
    // (For per-sequence updates, consider merging the update with the specific placeholder later.)
    if (key === "global" && statusData.percentage === 100 && statusData.result) {
      setFinalResults(statusData.result);
      sessionStorage.setItem("results", JSON.stringify(statusData.result));
    }
  }, []);

  // Initialize the SSE connection using the jobId stored in sessionStorage.
  useEffect(() => {
    const jobId = sessionStorage.getItem("jobId");
    if (jobId) {
      const eventSource = monitorProtocolProgress(jobId, onStatusUpdate);
      return () => eventSource.close();
    }
  }, [onStatusUpdate]);

  return (
    <div className="output-container">
      {progress.global ? (
        <div className="global-progress">
          <p>
            Global Progress: {progress.global.percentage}% — {progress.global.message}
          </p>
        </div>
      ) : (
        // If no SSE update has come in yet, display the initial message.
        initialMessage && (
          <div className="global-progress">
            <p>{initialMessage}</p>
          </div>
        )
      )}
      {dataToDisplay && dataToDisplay.length > 0 ? (
        <Results data={dataToDisplay} progress={progress} />
      ) : (
        <p className="initialization-message">Loading...</p>
      )}
    </div>
  );
}

export default ResultsPage;
