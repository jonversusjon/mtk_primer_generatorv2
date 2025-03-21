// FormPage.jsx
import { API_BASE_URL } from "../config/config.js";
import React, { useState, useEffect, useRef } from "react";
import { useNavigate } from "react-router-dom";
import Form from "../components/Form/Form";
import Sidebar from "../components/Form/Sidebar";
import { initiateProtocol } from "../api/api";
import useValidateForm from "../hooks/useValidateForm";
import { useFormUpdater } from "../hooks/useFormUpdater";
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
    kozak: "MTK",
    maxMutationsPerSite: 1,
    verboseMode: false,
    maxResults: "one",
  });

  const [loading, setLoading] = useState(true);
  const [configLoaded, setConfigLoaded] = useState(false);
  const [speciesLoaded, setSpeciesLoaded] = useState(false);
  const defaultsLoaded = configLoaded && speciesLoaded;

  const [processing, setProcessing] = useState(false);
  const [progressStatus, setProgressStatus] = useState({
    message: "",
    percentage: 0,
    step: "",
  });
  const [error, setError] = useState(null);
  const settingsToggleRef = useRef(null);
  const [activeTabIndex, setActiveTabIndex] = useState(0);
  const navigate = useNavigate();
  // Note: eventSourceRef is no longer needed here as SSE is managed in ResultsPage.

  const { updateSettings, updateFormInput } = useFormUpdater(setFormData);

  // Log formData changes for debugging.
  useEffect(() => {
    console.log("FormPage formData:", formData);
  }, [formData]);

  // Fetch initial configuration.
  useEffect(() => {
    const fetchConfig = async () => {
      try {
        const savedData = sessionStorage.getItem("formData");
        if (savedData) {
          const parsedData = JSON.parse(savedData);
          console.log("Loaded formData from sessionStorage:", parsedData);
          // Normalize sequences if needed.
          if (parsedData.sequencesToDomesticate) {
            parsedData.sequencesToDomesticate = parsedData.sequencesToDomesticate.map((seq) => ({
              ...seq,
              sequence: Array.isArray(seq.sequence) ? seq.sequence.join("") : seq.sequence,
            }));
          }
          setFormData((prev) => ({ ...prev, ...parsedData }));
        } else {
          const response = await fetch(`${API_BASE_URL}/config`);
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
                    sequence: Array.isArray(seq.sequence) ? seq.sequence.join("") : seq.sequence,
                  }))
                : [defaultSequence],
              species: data.species || "",
              kozak: data.kozak || "",
              maxMutationsPerSite: data.maxMutationsPerSite !== undefined ? data.maxMutationsPerSite : null,
              maxResults: data.maxResults !== undefined ? data.maxResults : null,
              verboseMode: data.verboseMode !== undefined ? data.verboseMode : null,
            };
          }
          console.log("New formData set from API:", newData);
          setFormData((prev) => ({ ...prev, ...newData }));
        }
      } catch (err) {
        console.error("Error fetching defaults from API:", err);
        setFormData((prev) => ({
          ...prev,
          sequencesToDomesticate: [defaultSequence],
        }));
      } finally {
        setLoading(false);
        setConfigLoaded(true);
      }
    };
    fetchConfig();
  }, []);

  // Fetch available species and update formData.
  useEffect(() => {
    const fetchSpecies = async () => {
      try {
        const response = await fetch(`${API_BASE_URL}/species`);
        const speciesData = await response.json();
        console.log("Fetched species from API:", speciesData);
        setFormData((prev) => ({
          ...prev,
          availableSpecies: speciesData.species,
          species: prev.species || (speciesData.species.length > 0 ? speciesData.species[0] : ""),
          kozak: prev.kozak || (speciesData.species.length > 0 ? speciesData.species[0] : ""),
          maxMutationsPerSite:
            prev.maxMutationsPerSite !== null && prev.maxMutationsPerSite !== undefined
              ? prev.maxMutationsPerSite
              : 1,
          maxResults:
            prev.maxResults !== null && prev.maxResults !== undefined ? prev.maxResults : "one",
          verboseMode:
            prev.verboseMode !== null && prev.verboseMode !== undefined ? prev.verboseMode : false,
        }));
      } catch (err) {
        console.error("Error fetching species:", err);
        setFormData((prev) => ({ ...prev, availableSpecies: [] }));
      } finally {
        setSpeciesLoaded(true);
      }
    };
    fetchSpecies();
  }, []);

  // Ensure species is set once both config and species have loaded.
  useEffect(() => {
    if (speciesLoaded && configLoaded && formData.availableSpecies.length > 0 && !formData.species) {
      setFormData((prev) => ({
        ...prev,
        species: prev.availableSpecies[0],
      }));
    }
  }, [speciesLoaded, configLoaded, formData.availableSpecies, formData.species]);

  // Validate the form.
  const { errors, isValid } = useValidateForm(formData, defaultsLoaded);
  const errorsBySequence = getErrorsBySequence(errors, formData.sequencesToDomesticate.length);

  // New handleFormSubmit: immediately initiate the protocol and set placeholders.
  const handleFormSubmit = async (data) => {
    console.log("Form submitted with data:", data);
    setProcessing(true);
    setError(null);
    setProgressStatus({
      message: "Initializing protocol generation...",
      percentage: 0,
      step: "init",
    });

    try {
      // 1. Store the submitted form data in sessionStorage.
      sessionStorage.setItem("formData", JSON.stringify(data));

      // 2. Generate a jobId if not provided.
      const jobId = data.jobId || Date.now().toString();
      const dataWithJobId = { ...formData, jobId };

      // 3. Save the jobId in sessionStorage for later use on the Results page.
      sessionStorage.setItem("jobId", jobId);

      // 4. Initiate the protocol to get the initial response.
      console.log("Calling initiateProtocol API helper");
      const initialData = await initiateProtocol(dataWithJobId);
      console.log("Initial protocol response:", initialData);

      // 5. Save the initial message from the server (or use a default).
      if (initialData.message) {
        sessionStorage.setItem("initialMessage", initialData.message);
      } else {
        sessionStorage.setItem("initialMessage", "Primer design started...");
      }

      // 6. Build placeholder objects from data.sequencesToDomesticate.
      const placeholders = data.sequencesToDomesticate.map((seq, idx) => ({
        id: idx,
        placeholder: true,
        sequence: seq.sequence,
        primerName: seq.primerName || `Sequence ${idx + 1}`,
      }));

      // 7. Store placeholders in sessionStorage and update top-level results.
      sessionStorage.setItem("results", JSON.stringify(placeholders));
      setResults(placeholders);

      // 8. Navigate immediately to the Results page.
      navigate("/results");
    } catch (err) {
      const errorMessage = err.message || "An error occurred while generating the protocol";
      console.error("Error in handleFormSubmit:", errorMessage, err);
      setError(errorMessage);
    } finally {
      setProcessing(false);
    }
  };

  if (loading || !defaultsLoaded) {
    return <div className="initialization-message">Getting things ready...</div>;
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
          updateSettings={updateSettings}
          formData={formData}
        />
        <div className="form-container" style={{ flex: 1, paddingLeft: "20px" }}>
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
            updateFields={updateFormInput}
            showSettings={showSettings}
            setShowSettings={setShowSettings}
            settingsToggleRef={settingsToggleRef}
            errors={errors}
            isValid={isValid}
            initialized={defaultsLoaded}
            activeTabIndex={activeTabIndex}
            setActiveTabIndex={setActiveTabIndex}
          />
        </div>
      </div>
    </div>
  );
}

export default FormPage;
