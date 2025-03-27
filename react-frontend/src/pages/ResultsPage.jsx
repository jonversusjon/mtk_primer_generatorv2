// ResultsPage.jsx
import React, { useEffect, useMemo } from "react";
import { useNavigate } from "react-router-dom";
import ResultTabs from "../components/Results/ResultTabs";

function ResultsPage({ results }) {
  const navigate = useNavigate();
  // Get the jobId from sessionStorage (or from any other source you use)
  const jobId = sessionStorage.getItem("jobId");

  // Build placeholder data from stored form data if needed.
  const placeholders = useMemo(() => {
    const savedFormData = sessionStorage.getItem("formData");
    if (savedFormData) {
      const parsed = JSON.parse(savedFormData);
      return (
        parsed.sequencesToDomesticate?.map((seq, i) => ({
          id: i, // used as sequenceIdx later
          placeholder: true,
          sequence: seq.sequence,
          primerName: seq.primerName || `Sequence ${i + 1}`,
        })) || []
      );
    }
    return [];
  }, []);

  const dataToDisplay = results || placeholders;

  // Redirect to form if no data is available.
  useEffect(() => {
    if (!dataToDisplay?.length) {
      console.log("No data found — redirecting to form");
      navigate("/");
    }
  }, [dataToDisplay, navigate]);

  return (
    <div className="output-container">
      {dataToDisplay?.length ? (
        // Pass jobId along with the data so that each ResultTab can receive both
        <ResultTabs results={dataToDisplay} jobId={jobId} />
      ) : (
        <p className="initialization-message">Loading...</p>
      )}
      <div className="h-16" />
    </div>
  );
}

export default ResultsPage;
