// ResultsPage.jsx
import React, { useEffect, useMemo } from "react";
import { useNavigate } from "react-router-dom";
import ResultTabs from "../components/Results//ResultTabs";

function ResultsPage({ results }) {
  const navigate = useNavigate();

  // Build placeholder data from stored form data if needed.
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
        <ResultTabs results={dataToDisplay} />
      ) : (
        <p className="initialization-message">Loading...</p>
      )}
      <div className="h-16" />
    </div>
  );
}

export default ResultsPage;
