// src/pages/FormPage.jsx
import React, { useState, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import Form from "../components/Form/Form";
import { generateProtocol } from "../api/api";
import { getDefaultValues } from "../config/defaultParameters"; // Hardcoded fallback

function FormPage({ showSettings, setShowSettings, setResults }) {
  const [formData, setFormData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const navigate = useNavigate();

  // default values cascade: sessionStorage > Flask API > hardcoded defaults
  useEffect(() => {
    // Step 1: Check session storage first
    const savedData = sessionStorage.getItem("formData");
    if (savedData) {
      setFormData(JSON.parse(savedData));
    } else {
      // Step 2: If no session data, fetch defaults from Flask
      setLoading(true);
      fetch(`http://localhost:5000/api/config/development`)
        .then((res) => res.json())
        .then((data) => {
          // Step 3: If API returns empty or fails, fall back to hardcoded defaults
          if (!data || Object.keys(data).length === 0) {
            console.warn("API returned empty config, using hardcoded defaults");
            setFormData(getDefaultValues());
          } else {
            setFormData(data);
          }
        })
        .catch((err) => {
          console.error("Error fetching defaults from Flask:", err);
          setFormData(getDefaultValues()); // Ensure fallback always works
        })
        .finally(() => setLoading(false));
    }
  }, []);

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
