import React, { useState, useEffect, useCallback, useMemo } from "react";
import { useNavigate } from "react-router-dom";
import Results from "../components/Results/Results";
import { monitorProtocolProgress } from "../api/api";

function ResultsPage({ results }) {
  const navigate = useNavigate();

  const [finalResults, setFinalResults] = useState(() => {
    const saved = sessionStorage.getItem("results");
    return saved ? JSON.parse(saved) : null;
  });

  const [initialMessage] = useState(() => {
    const msg =
      sessionStorage.getItem("initialMessage") || "Primer design started...";
    console.log("Initial message:", msg);
    return msg;
  });

  const placeholders = useMemo(() => {
    const savedFormData = sessionStorage.getItem("formData");
    if (savedFormData) {
      const parsed = JSON.parse(savedFormData);
      return (
        parsed.sequencesToDomesticate?.map((seq, i) => ({
          id: i,
          placeholder: true,
          sequence: seq.sequence,
          primerName: seq.primerName || `Sequence ${i + 1}`,
        })) || []
      );
    }
    return [];
  }, []);

  const [progress, setProgress] = useState({});

  const dataToDisplay = finalResults || results || placeholders;

  useEffect(() => {
    if (!dataToDisplay?.length) {
      console.log("No data found — redirecting to form");
      navigate("/");
    }
  }, [dataToDisplay, navigate]);

  const onStatusUpdate = useCallback((statusData) => {
    const key = statusData.sequenceId ?? "global";
    console.log("Received SSE update:", statusData);

    setProgress((prev) => ({ ...prev, [key]: statusData }));

    if (
      key === "global" &&
      statusData.percentage === 100 &&
      statusData.result
    ) {
      console.log("Global complete — saving final results:", statusData.result);
      setFinalResults(statusData.result);
      sessionStorage.setItem("results", JSON.stringify(statusData.result));
    }
  }, []);

  useEffect(() => {
    const jobId = sessionStorage.getItem("jobId");
    if (jobId) {
      console.log("Initializing SSE for jobId:", jobId);
      const eventSource = monitorProtocolProgress(jobId, onStatusUpdate);
      return () => {
        console.log("Closing SSE connection");
        eventSource.close();
      };
    } else {
      console.warn("No jobId in sessionStorage — cannot monitor progress");
    }
  }, [onStatusUpdate]);

  return (
    <div className="output-container">
      {progress.global ? (
        <div className="global-progress">
          <p>
            Global Progress: {progress.global.percentage}% —{" "}
            {progress.global.message}
          </p>
        </div>
      ) : (
        initialMessage && (
          <div className="global-progress">
            <p>{initialMessage}</p>
          </div>
        )
      )}
      {dataToDisplay?.length ? (
        <Results data={dataToDisplay} progress={progress} />
      ) : (
        <p className="initialization-message">Loading...</p>
      )}
    </div>
  );
}

export default ResultsPage;
