// ResultsPage.jsx
import React, { useEffect, useMemo } from "react";
import { useNavigate } from "react-router-dom";
import Results from "../components/Results/Results";

function ResultsPage({ results }) {
  const navigate = useNavigate();

  // const [initialMessage] = useState(() => {
  //   const msg =
  //     sessionStorage.getItem("initialMessage") || "Primer design started...";
  //   console.log("Initial message:", msg);
  //   return msg;
  // });

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

  useEffect(() => {
    if (!dataToDisplay?.length) {
      console.log("No data found — redirecting to form");
      navigate("/");
    }
  }, [dataToDisplay, navigate]);

  return (
    <div className="output-container">
      {dataToDisplay?.length ? (
        <Results data={dataToDisplay} />
      ) : (
        <p className="initialization-message">Loading...</p>
      )}
      <div className="h-16" />
    </div>
  );
}

export default ResultsPage;
