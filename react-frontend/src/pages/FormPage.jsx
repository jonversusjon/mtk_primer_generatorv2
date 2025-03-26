import { API_BASE_URL } from "../config/config.js";
import React, { useState, useEffect, useRef } from "react";
import { useNavigate } from "react-router-dom";
import Form from "../components/Form/Form";
import Sidebar from "../components/Form/Sidebar";
import { submitProtocol } from "../api/api";
import useValidateForm from "../hooks/useValidateForm";
import { useFormUpdater } from "../hooks/useFormUpdater";
import "../styles/Form.css";

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

  const [prefillLoaded, setPrefillLoaded] = useState(false);
  const [speciesLoaded, setSpeciesLoaded] = useState(false);
  const [processing, setProcessing] = useState(false);
  const [error, setError] = useState(null);
  const settingsToggleRef = useRef(null);
  const [activeTabIndex, setActiveTabIndex] = useState(0);
  const navigate = useNavigate();
  const { updateSettings, updateFormInput } = useFormUpdater(setFormData);

  // Load persisted → dummy → default
  useEffect(() => {
    const initializeFormData = async () => {
      try {
        const saved = sessionStorage.getItem("formData");
        if (saved) {
          setFormData(JSON.parse(saved));
        } else {
          // 1️⃣ Fetch dummy data
          const dummyResp = await fetch(`${API_BASE_URL}/dummy`);
          const dummy = dummyResp.ok ? await dummyResp.json() : {};
  
          // 2️⃣ Fetch species list
          const speciesResp = await fetch(`${API_BASE_URL}/species`);
          const { species: speciesList = [] } = await speciesResp.json();
  
          // 3️⃣ Merge defaults
          setFormData({
            ...dummy,
            availableSpecies: speciesList,
            species: dummy.species || speciesList[0] || "",
            kozak: dummy.kozak || speciesList[0] || "",
          });
        }
      } catch (err) {
        console.error("Error loading prefill data", err);
      } finally {
        setPrefillLoaded(true);
        setSpeciesLoaded(true);
      }
    };
    initializeFormData();
  }, []);
  

  // Load species dropdown options
  useEffect(() => {
    const fetchSpecies = async () => {
      try {
        const resp = await fetch(`${API_BASE_URL}/species`);
        const { species } = await resp.json();
        setFormData((prev) => ({
          ...prev,
          availableSpecies: species,
          species: prev.species || species[0] || "",
          kozak: prev.kozak || species[0] || "",
        }));
      } catch {
        console.error("Error fetching species");
      } finally {
        setSpeciesLoaded(true);
      }
    };
    fetchSpecies();
  }, []);

  // Validation
  const { errors, isValid } = useValidateForm(
    formData,
    prefillLoaded && speciesLoaded
  );

  const handleFormSubmit = async (data) => {
    setProcessing(true);
    setError(null);

    try {
      sessionStorage.setItem("formData", JSON.stringify(data));
      const { jobId } = await submitProtocol(data);
      sessionStorage.setItem("jobId", jobId);

      const placeholders = data.sequencesToDomesticate.map((seq, idx) => ({
        id: idx,
        placeholder: true,
        sequence: seq.sequence,
        primerName: seq.primerName || `Sequence ${idx + 1}`,
      }));
      sessionStorage.setItem("results", JSON.stringify(placeholders));
      sessionStorage.setItem("jobId", jobId);
      setResults(placeholders);
      navigate("/results");
    } catch (err) {
      console.error("Submit error:", err);
      setError(err.message || "Failed to start protocol generation");
    } finally {
      setProcessing(false);
    }
  };

  if (!prefillLoaded || !speciesLoaded) {
    return <div className="initialization-message">Loading form…</div>;
  }

  return (
    <div className="form-page-container">
      <header className="form-header">
        <h2>Primer Design Form</h2>
      </header>
      <div className="form-layout">
        <Sidebar
          sequences={formData.sequencesToDomesticate}
          errorsBySequence={errors.sequencesToDomesticate || []}
          onSelectTab={setActiveTabIndex}
          activeTabIndex={activeTabIndex}
          showSettings={showSettings}
          setShowSettings={setShowSettings}
          settingsToggleRef={settingsToggleRef}
          updateSettings={updateSettings}
          formData={formData}
        />
        <main className="form-container">
          {error && <div className="alert alert-danger">{error}</div>}
          {processing && <p>Processing…</p>}
          <Form
            onSubmit={handleFormSubmit}
            formData={formData}
            updateFields={updateFormInput}
            showSettings={showSettings}
            initialized={prefillLoaded && speciesLoaded}
            errorsBySequence={errors.sequencesToDomesticate || []}
            isValid={isValid}
            activeTabIndex={activeTabIndex}
            setActiveTabIndex={setActiveTabIndex}
          />
        </main>
      </div>
      <div className="h-16" />
    </div>
  );
}

export default FormPage;
