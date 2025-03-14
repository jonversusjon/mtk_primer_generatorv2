import React, { useState, useEffect, useRef } from "react";
import { useNavigate } from "react-router-dom";
import Form from "../components/Form/Form";
import Sidebar from "../components/Form/Sidebar";
import { generateProtocol } from "../api/api";
import useValidateForm from "../hooks/useValidateForm";
import "../styles/Form.css";

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
  const [initialized, setInitialized] = useState(false); // NEW: indicates data is fully loaded
  const [error, setError] = useState(null);
  const settingsToggleRef = useRef(null);
  const [activeTabIndex, setActiveTabIndex] = useState(0);
  const navigate = useNavigate();

  useEffect(() => {
    console.log("FormPage formData:", formData);
  }, [formData]);

  useEffect(() => {
    const fetchConfig = async () => {
      try {
        const savedData = sessionStorage.getItem("formData");
        if (savedData) {
          const parsedData = JSON.parse(savedData);
          console.log("Loaded formData from sessionStorage:", parsedData); // DEBUG log
          // Normalize sequences if needed
          if (parsedData.sequencesToDomesticate) {
            parsedData.sequencesToDomesticate =
              parsedData.sequencesToDomesticate.map((seq) => ({
                ...seq,
                sequence: Array.isArray(seq.sequence)
                  ? seq.sequence.join("")
                  : seq.sequence,
              }));
          }
          setFormData(parsedData);
        } else {
          const response = await fetch("http://localhost:5000/api/config");
          const data = await response.json();
          console.log("Fetched formData from API:", data); // DEBUG log
          let newData;
          if (!data || Object.keys(data).length === 0) {
            console.warn("API returned empty config! Using fallback defaults.");
            newData = { sequencesToDomesticate: [defaultSequence] };
          } else {
            newData = {
              ...data,
              sequencesToDomesticate: data.sequencesToDomesticate
                ? data.sequencesToDomesticate.map((seq) => ({
                    ...seq,
                    sequence: Array.isArray(seq.sequence)
                      ? seq.sequence.join("")
                      : seq.sequence,
                  }))
                : [defaultSequence],
            };
          }
          console.log("New formData set from API:", newData); // DEBUG log
          setFormData(newData);
        }
        setInitialized(true);
      } catch (err) {
        console.error("Error fetching defaults from Flask:", err);
        setFormData({ sequencesToDomesticate: [defaultSequence] });
        setInitialized(true);
      } finally {
        setLoading(false);
      }
    };
    fetchConfig();
  }, []);

  // Now we pass `initialized` so that validation doesn't run until formData is ready.
  const { errors, isValid } = useValidateForm(formData, initialized);
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
      <div className="form-header" style={{ width: "100%" }}>
        <h2 className="primer-form-title">Primer Design Form</h2>
      </div>
      <div style={{ display: "flex" }}>
        <Sidebar
          sequences={formData.sequencesToDomesticate}
          errorsBySequence={errorsBySequence}
          onSelectTab={setActiveTabIndex}
          activeTabIndex={activeTabIndex}
          settingsToggleRef={settingsToggleRef}
          setShowSettings={setShowSettings}
          showSettings={showSettings}
          formData={formData}
          setFormData={setFormData}
        />
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
