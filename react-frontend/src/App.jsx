import React, { useState } from "react";
import AppHeader from "./components/AppHeader";
import Form from "./components/Form/Form";
import Results from "./components/Results/Results";
import { useDarkMode } from "./hooks/useDarkMode";
import { generateProtocol } from "./api/api";
import "./styles/app.css";

function App() {
  const [darkMode, toggleDarkMode] = useDarkMode();
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const handleFormSubmit = async (formData) => {
    setLoading(true);
    setError(null);

    try {
      // This will be implemented in the API section
      const response = await generateProtocol(formData);
      setResults(response);
    } catch (err) {
      setError(
        err.message || "An error occurred while generating the protocol"
      );
      console.error("Error generating protocol:", err);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className={`app ${darkMode ? "dark-mode" : ""}`}>
      <AppHeader darkMode={darkMode} toggleDarkMode={toggleDarkMode} />

      <div className="app-container">
        <div className="form-container">
          <Form onSubmit={handleFormSubmit} isSubmitting={loading} />
        </div>

        <div className="output-container">
          {error && (
            <div id="error-container" className="alert alert-danger">
              {error}
            </div>
          )}

          {loading ? (
            <div className="loading-indicator">
              <div className="spinner-container">
                <div className="spinner-border" role="status">
                  <span className="visually-hidden">Processing...</span>
                </div>
                <span className="loading-text">Processing your request...</span>
              </div>
            </div>
          ) : (
            results && <Results data={results} />
          )}
        </div>
      </div>
    </div>
  );
}

export default App;
