import React, { useState, useEffect, useCallback, useMemo } from "react";
import { useNavigate } from "react-router-dom";
import Results from "../components/Results/Results";
import { monitorProtocolProgress } from "../api/api";

// Custom hook to manage SSE connection
const useJobProgress = (onStatusUpdate) => {
  useEffect(() => {
    const jobId = sessionStorage.getItem("jobId");
    if (jobId) {
      const eventSource = monitorProtocolProgress(jobId, onStatusUpdate);
      return () => eventSource.close();
    }
  }, [onStatusUpdate]);
};

function ResultsPage({ results }) {
  const navigate = useNavigate();

  // Load final results from sessionStorage (if any) into state.
  const [finalResults, setFinalResults] = useState(() => {
    const savedResults = sessionStorage.getItem("results");
    return savedResults ? JSON.parse(savedResults) : null;
  });

  // Load formData from sessionStorage and create placeholder objects.
  // This creates one placeholder per sequence submitted.
  const placeholders = useMemo(() => {
    const savedFormData = sessionStorage.getItem("formData");
    console.log("Saved form data:", savedFormData);
    if (savedFormData) {
      const parsedFormData = JSON.parse(savedFormData);
      if (parsedFormData.sequencesToDomesticate) {
        return parsedFormData.sequencesToDomesticate.map((seq, index) => ({
          // You can adjust the placeholder object structure as needed.
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
  const [progress, setProgress] = useState({});

  // Use either the provided results prop, finalResults from SSE, or placeholders.
  const dataToDisplay = finalResults || results || placeholders;

  // If no results or form data, redirect to form.
  useEffect(() => {
    if (!dataToDisplay || dataToDisplay.length === 0) {
      navigate("/");
    }
  }, [dataToDisplay, navigate]);

  // Status update callback from SSE events.
  const onStatusUpdate = useCallback((statusData) => {
    // Update progress state based on a key, which might be sequence-specific or global.
    const key =
      statusData.sequenceId !== undefined ? statusData.sequenceId : "global";
    setProgress((prev) => ({ ...prev, [key]: statusData }));

    // If final data arrives (100% progress), store it in local state.
    if (statusData.percentage === 100 && statusData.result) {
      console.log("Final protocol data received via SSE:", statusData.result);
      setFinalResults(statusData.result);
      // Optionally, you can also save the result to sessionStorage.
      sessionStorage.setItem("results", JSON.stringify(statusData.result));
    }
  }, []);

  // Initialize the SSE connection if a jobId exists.
  useJobProgress(onStatusUpdate);

  return (
    <div className="output-container">
      {progress.global && (
        <div className="global-progress">
          <p>
            Global Progress: {progress.global.percentage}% —{" "}
            {progress.global.message}
          </p>
        </div>
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
