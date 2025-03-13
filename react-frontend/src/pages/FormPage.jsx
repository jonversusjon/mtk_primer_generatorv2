import React, { useState, useEffect, useRef } from "react";
import { useNavigate } from "react-router-dom";
import Form from "../components/Form/Form";
import Sidebar from "../components/Form/Sidebar";
import { generateProtocol } from "../api/api";
import useValidateForm from "../hooks/useValidateForm";
import "../styles/Form.css";

// Helper function to group validation errors by sequence index
const getErrorsBySequence = (errors, count) => {
  const errorsBySequence = Array.from({ length: count }, () => []);
  Object.entries(errors).forEach(([key, message]) => {
    const match = key.match(/sequencesToDomesticate\[(\d+)\]/);
    if (match) {
      const index = parseInt(match[1], 10);
      if (index < count) {
        errorsBySequence[index].push(message);
      }
    }
  });
  return errorsBySequence;
};

const defaultSequence = {
  sequence: "",
  primerName: "",
  mtkPartLeft: "",
  mtkPartRight: "",
};

function FormPage({ showSettings, setShowSettings, setResults }) {
  const [formData, setFormData] = useState({
    sequencesToDomesticate: [defaultSequence],
  });
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  // Ref for the settings button
  const settingsToggleRef = useRef(null);

  // Manage which tab is active
  const [activeTabIndex, setActiveTabIndex] = useState(0);

  const navigate = useNavigate();

  useEffect(() => {
    const fetchConfig = async () => {
      try {
        const savedData = sessionStorage.getItem("formData");
        if (savedData) {
          setFormData(JSON.parse(savedData));
        } else {
          const response = await fetch("http://localhost:5000/api/config");
          const data = await response.json();
          if (!data || Object.keys(data).length === 0) {
            console.warn("API returned empty config! Using fallback defaults.");
            setFormData({ sequencesToDomesticate: [defaultSequence] });
          } else {
            setFormData({
              ...data,
              sequencesToDomesticate: data.sequencesToDomesticate || [
                defaultSequence,
              ],
            });
          }
        }
      } catch (err) {
        console.error("Error fetching defaults from Flask:", err);
        setFormData({ sequencesToDomesticate: [defaultSequence] });
      } finally {
        setLoading(false);
      }
    };
    fetchConfig();
  }, []);

  // Run validation in the parent
  const { errors, isValid } = useValidateForm(formData);
  const errorsBySequence = getErrorsBySequence(
    errors,
    formData.sequencesToDomesticate.length
  );

  const handleFormSubmit = async (data) => {
    setLoading(true);
    setError(null);

    try {
      sessionStorage.setItem("formData", JSON.stringify(data));
      const response = await generateProtocol(data);
      sessionStorage.setItem("results", JSON.stringify(response));
      setResults(response);
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

  if (loading) {
    return <p>Loading defaults...</p>;
  }

  return (
    <div className="form-page-container">
      {/* Header that spans across the entire width */}
      <div className="form-header" style={{ width: "100%" }}>
        <h2 className="primer-form-title">Primer Design Form</h2>
      </div>

      {/* Content container for sidebar and form side by side */}
      <div style={{ display: "flex" }}>
        {/* Sidebar on the left */}
        <Sidebar
          sequences={formData.sequencesToDomesticate}
          errorsBySequence={errorsBySequence}
          onSelectTab={setActiveTabIndex}
          activeTabIndex={activeTabIndex}
        />

        {/* Form on the right */}
        <div
          className="form-container"
          style={{ flex: 1, paddingLeft: "20px" }}
        >
          {error && <div className="alert alert-danger">{error}</div>}

          <Form
            onSubmit={handleFormSubmit}
            formData={formData}
            setFormData={setFormData}
            showSettings={showSettings}
            setShowSettings={setShowSettings}
            settingsToggleRef={settingsToggleRef}
            errors={errors}
            isValid={isValid}
            activeTabIndex={activeTabIndex}
            setActiveTabIndex={setActiveTabIndex}
          />
        </div>
      </div>
    </div>
  );
}

export default FormPage;
