// FormPage.jsx
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
    availableSpecies: [],
    species: "",
  });
  const [loading, setLoading] = useState(true);
  const [processing, setProcessing] = useState(false);
  const [progressStatus, setProgressStatus] = useState({
    message: "",
    percentage: 0,
    step: "",
  });
  const [initialized, setInitialized] = useState(false);
  const [error, setError] = useState(null);
  const settingsToggleRef = useRef(null);
  const [activeTabIndex, setActiveTabIndex] = useState(0);
  const navigate = useNavigate();
  const eventSourceRef = useRef(null);

  useEffect(() => {
    console.log("FormPage formData:", formData);
  }, [formData]);

  // Fetch initial configuration
  useEffect(() => {
    const fetchConfig = async () => {
      try {
        const savedData = sessionStorage.getItem("formData");
        if (savedData) {
          const parsedData = JSON.parse(savedData);
          console.log("Loaded formData from sessionStorage:", parsedData);
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
          setFormData((prev) => ({ ...prev, ...parsedData }));
        } else {
          const response = await fetch("http://localhost:5000/api/config");
          const data = await response.json();
          console.log("Fetched formData from API:", data);
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
          console.log("New formData set from API:", newData);
          setFormData((prev) => ({ ...prev, ...newData }));
        }
        setInitialized(true);
      } catch (err) {
        console.error("Error fetching defaults from API:", err);
        setFormData((prev) => ({
          ...prev,
          sequencesToDomesticate: [defaultSequence],
        }));
        setInitialized(true);
      } finally {
        setLoading(false);
      }
    };
    fetchConfig();
  }, []);

  // Fetch available species and update formData
  useEffect(() => {
    const fetchSpecies = async () => {
      try {
        const response = await fetch("http://localhost:5000/api/species");
        const speciesData = await response.json();
        console.log("Fetched species from API:", speciesData);
        setFormData((prev) => ({
          ...prev,
          availableSpecies: speciesData.species,
          // If no species is selected, default to the first available species.
          species:
            prev.species ||
            (speciesData.species.length > 0 ? speciesData.species[0] : ""),
        }));
      } catch (err) {
        console.error("Error fetching species:", err);
        setFormData((prev) => ({ ...prev, availableSpecies: [] }));
      }
    };
    fetchSpecies();
  }, []);

  // Clean up event source when component unmounts
  useEffect(() => {
    return () => {
      if (eventSourceRef.current) {
        eventSourceRef.current.close();
      }
    };
  }, []);

  // Wait for the form to be initialized before running validations
  const { errors, isValid } = useValidateForm(formData, initialized);
  const errorsBySequence = getErrorsBySequence(
    errors,
    formData.sequencesToDomesticate.length
  );

  const handleFormSubmit = async (data) => {
    setProcessing(true);
    setError(null);
    setProgressStatus({
      message: "Initializing protocol generation...",
      percentage: 0,
      step: "init",
    });

    try {
      sessionStorage.setItem("formData", JSON.stringify(data));
      // Generate a job ID and include it in your data payload
      const jobId = Date.now().toString();

      // Call generateProtocol with the onStatusUpdate callback.
      // The API function will handle creating and managing the SSE connection.
      const result = await generateProtocol({ ...data, jobId }, (status) => {
        setProgressStatus(status);
      });

      // Optionally store the eventSource for cleanup if needed.
      if (result.eventSource) {
        eventSourceRef.current = result.eventSource;
      }

      sessionStorage.setItem("results", JSON.stringify(result));
      setResults(result);
      navigate("/results");
    } catch (err) {
      setError(
        err.message || "An error occurred while generating the protocol"
      );
      console.error("Error generating protocol:", err);
    } finally {
      setProcessing(false);
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

          {processing && (
            <div className="processing-status">
              <div className="progress">
                <div
                  className="progress-bar"
                  role="progressbar"
                  style={{ width: `${progressStatus.percentage}%` }}
                  aria-valuenow={progressStatus.percentage}
                  aria-valuemin="0"
                  aria-valuemax="100"
                >
                  {progressStatus.percentage}%
                </div>
              </div>
              <p>{progressStatus.message}</p>
              <p>Current step: {progressStatus.step}</p>
            </div>
          )}

          <Form
            onSubmit={handleFormSubmit}
            formData={formData}
            setFormData={setFormData}
            showSettings={showSettings}
            setShowSettings={setShowSettings}
            settingsToggleRef={settingsToggleRef}
            errors={errors}
            isValid={isValid}
            initialized={initialized}
            activeTabIndex={activeTabIndex}
            setActiveTabIndex={setActiveTabIndex}
          />
        </div>
      </div>
    </div>
  );
}

export default FormPage;
