// src/pages/ResultsPage.jsx
import React, { useState, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import Results from "../components/Results/Results";

// TODO: the returned data needs to be reformatted to match the expected format of the Results component
// TODO: initial return should check how many sequencesToDomesticate to make just one tabe per sequence

function ResultsPage({ results }) {
  const [localResults, setLocalResults] = useState(results);
  const navigate = useNavigate();

  useEffect(() => {
    // If no results in parent state, try sessionStorage
    if (!results) {
      const storedResults = sessionStorage.getItem("results");
      if (storedResults) {
        setLocalResults(JSON.parse(storedResults));
      } else {
        // If we truly have nothing, send user back to the form
        navigate("/");
      }
    }
  }, [results, navigate]);

  return (
    <div className="output-container">
      {localResults ? <Results data={localResults} /> : <p>Loading...</p>}
    </div>
  );
}

export default ResultsPage;
