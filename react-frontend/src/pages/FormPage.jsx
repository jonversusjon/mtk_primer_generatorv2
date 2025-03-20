// FormPage.jsx
import { API_BASE_URL } from "../config/config.js";

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
    kozak: "MTK",
    maxMutationsPerSite: 1,
    maxResults: 1,
    verboseMode: false,
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
  const eventSourceRef = useRef(null);

  // Log formData changes
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
                    sequence: Array.isArray(seq.sequence)
                      ? seq.sequence.join("")
                      : seq.sequence,
                  }))
                : [defaultSequence],
              species: data.species || "",
              kozak: data.kozak || "",
              maxMutationsPerSite:
                data.maxMutationsPerSite !== undefined
                  ? data.maxMutationsPerSite
                  : null,
              max_results:
                data.max_results !== undefined ? data.max_results : null,
              verboseMode:
                data.verboseMode !== undefined ? data.verboseMode : null,
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

  // Fetch available species and update formData
  useEffect(() => {
    const fetchSpecies = async () => {
      try {
        const response = await fetch(`${API_BASE_URL}/species`);
        const speciesData = await response.json();
        console.log("Fetched species from API:", speciesData);
        setFormData((prev) => ({
          ...prev,
          availableSpecies: speciesData.species,
          species:
            prev.species ||
            (speciesData.species.length > 0 ? speciesData.species[0] : ""),
          kozak:
            prev.kozak ||
            (speciesData.species.length > 0 ? speciesData.species[0] : ""),
          maxMutationsPerSite:
            prev.maxMutationsPerSite !== null &&
            prev.maxMutationsPerSite !== undefined
              ? prev.maxMutationsPerSite
              : 1,
          max_results:
            prev.max_results !== null && prev.max_results !== undefined
              ? prev.max_results
              : 1,
          verboseMode:
            prev.verboseMode !== null && prev.verboseMode !== undefined
              ? prev.verboseMode
              : false,
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
    if (
      speciesLoaded &&
      configLoaded &&
      formData.availableSpecies.length > 0 &&
      !formData.species
    ) {
      setFormData((prev) => ({
        ...prev,
        species: prev.availableSpecies[0],
      }));
    }
  }, [
    speciesLoaded,
    configLoaded,
    formData.availableSpecies,
    formData.species,
  ]);

  // Clean up event source on component unmount.
  useEffect(() => {
    return () => {
      if (eventSourceRef.current) {
        console.log("Closing event source on unmount");
        eventSourceRef.current.close();
      }
    };
  }, []);

  // Manual implementation of protocol generation for better logging
  const manualGenerateProtocol = async (data) => {
    console.log("Starting manual protocol generation with data:", data);

    try {
      // Step 1: Call generate_protocol endpoint
      const response = await fetch(`${API_BASE_URL}/generate_protocol`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(data),
      });

      if (!response.ok) {
        const errorText = await response.text();
        console.error("Error response from generate_protocol:", errorText);
        throw new Error(`API error: ${response.status} ${errorText}`);
      }

      const initialData = await response.json();
      console.log("Initial protocol response:", initialData);

      const jobId = initialData.jobId;

      // Step 2: Set up SSE connection for status updates
      const eventSource = new EventSource(`${API_BASE_URL}/status/${jobId}`);
      eventSourceRef.current = eventSource;

      return new Promise((resolve, reject) => {
        eventSource.onmessage = (event) => {
          try {
            const status = JSON.parse(event.data);
            console.log("SSE status update:", status);

            // Update progress in UI
            setProgressStatus(status);

            // Check if job is complete
            if (status.percentage === 100) {
              console.log("Job complete! Results:", status.result);
              eventSource.close();
              resolve(status.result);
            }

            // Check if job failed
            if (status.percentage === -1) {
              const errorMsg = status.message || "Job failed";
              console.error("Job failed:", errorMsg);
              eventSource.close();
              reject(new Error(errorMsg));
            }
          } catch (err) {
            console.error("Error parsing SSE message:", err, event.data);
            reject(err);
          }
        };

        eventSource.onerror = (err) => {
          console.error("SSE connection error:", err);
          eventSource.close();
          reject(new Error("Connection to status stream failed"));
        };
      });
    } catch (err) {
      console.error("Error in manual protocol generation:", err);
      throw err;
    }
  };

  // Wait for the form to be initialized before running validations
  const { errors, isValid } = useValidateForm(formData, defaultsLoaded);
  const errorsBySequence = getErrorsBySequence(
    errors,
    formData.sequencesToDomesticate.length
  );

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
      // Store the entire form data in sessionStorage.
      sessionStorage.setItem("formData", JSON.stringify(data));

      // Generate a jobId if not already provided.
      const jobId = data.jobId || Date.now().toString();
      // Add jobId to the form data.
      const dataWithJobId = { ...data, jobId };

      // Store the jobId in sessionStorage for later use on the Results page.
      sessionStorage.setItem("jobId", jobId);

      // Use the existing generateProtocol function (with fallback to manualGenerateProtocol)
      let result;
      try {
        console.log("Using generateProtocol API helper");
        result = await generateProtocol(dataWithJobId, (status) => {
          console.log("Status update from API helper:", status);
          setProgressStatus(status);
        });
      } catch (apiError) {
        console.warn(
          "generateProtocol API helper failed, falling back to manual implementation:",
          apiError
        );
        result = await manualGenerateProtocol(dataWithJobId);
      }

      console.log("Protocol generation successful, final result:", result);

      // Save protocol results in sessionStorage.
      sessionStorage.setItem("results", JSON.stringify(result));
      setResults(result);
      // Navigate to the /results page after job has started.
      navigate("/results");
    } catch (err) {
      const errorMessage =
        err.message || "An error occurred while generating the protocol";
      console.error("Error in handleFormSubmit:", errorMessage, err);
      setError(errorMessage);
    } finally {
      setProcessing(false);
    }
  };

  if (loading || !defaultsLoaded) {
    return (
      <div className="initialization-message">Getting things ready...</div>
    );
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
