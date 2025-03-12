// src/pages/FormPage.jsx
import React, { useState } from "react";
import { useNavigate } from "react-router-dom";
import Form from "../components/Form/Form";
import { generateProtocol } from "../api/api";

function FormPage({ showSettings, setShowSettings, setResults }) {
  const [formData, setFormData] = useState(() => {
    // Recover previously typed form data from sessionStorage
    const savedData = sessionStorage.getItem("formData");
    return savedData ? JSON.parse(savedData) : {};
  });
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const navigate = useNavigate();

  const handleFormSubmit = async (data) => {
    setLoading(true);
    setError(null);

    try {
      // Persist the form data so user doesn't lose it on refresh
      sessionStorage.setItem("formData", JSON.stringify(data));

      // Call your Flask API
      const response = await generateProtocol(data);
      console.log("📢 Received API Response from Flask:", response);

      // Store results in sessionStorage so it can be recovered later
      sessionStorage.setItem("results", JSON.stringify(response));

      // Also store in parent state for immediate usage
      setResults(response);

      // Navigate to results page
      navigate("/results");
    } catch (err) {
      setError(
        err.message || "An error occurred while generating the protocol"
      );
      console.error("Error generating protocol:", err);
    } finally {
      setLoading(false);
    }
  };

  // If you want to track form field changes on the fly:
  const handleFormChange = (updatedData) => {
    setFormData(updatedData);
  };

  return (
    <div className="form-container">
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
        <Form
          onSubmit={handleFormSubmit}
          onChange={handleFormChange}
          formData={formData}
          showSettings={showSettings}
          setShowSettings={setShowSettings}
        />
      )}
    </div>
  );
}

export default FormPage;
