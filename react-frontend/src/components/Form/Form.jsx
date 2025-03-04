// Form.js
import React, { useState, useEffect } from "react";
import TemplateSequence from "./TemplateSequence";
import SequenceTabs from "./SequenceTabs";
import Settings from "./Settings";
import { fetchAvailableSpecies } from "../../api/api";
import { defaultParameters } from "../../config/defaultParameters";
import useValidateForm from "../../hooks/useValidateForm";
import "../../styles/form.css";

function getDefaultFormData() {
  console.log("Loading defaults with parameters:", defaultParameters);
  const initialSequences = [];

  if (
    defaultParameters.sequencesToDomesticate &&
    defaultParameters.sequencesToDomesticate.length > 0
  ) {
    defaultParameters.sequencesToDomesticate.forEach((seq, index) => {
      initialSequences.push({
        sequence: seq,
        primerName: defaultParameters.primerNames[index] || "",
        mtkPart: defaultParameters.mtkPartNums[index] || "",
      });
    });
  } else {
    // Default to one empty sequence
    initialSequences.push({ sequence: "", primerName: "", mtkPart: "" });
  }

  return {
    templateSequence: defaultParameters.templateSequence || "",
    numSequences: initialSequences.length,
    species: "", // We'll set this later once species is fetched
    kozak: "MTK",
    max_mut_per_site: 1,
    verbose_mode: true,
    sequences: initialSequences,
  };
}

function Form({ onSubmit, isSubmitting }) {
  // State for species loading and errors
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [species, setSpecies] = useState([]);

  // Initialize form state with defaults
  const [formData, setFormData] = useState(getDefaultFormData());

  // Settings panel state
  const [showSettings, setShowSettings] = useState(false);

  // Load available species when the component mounts
  useEffect(() => {
    const loadSpecies = async () => {
      try {
        setLoading(true);
        setError(null);
        const speciesData = await fetchAvailableSpecies();
        setSpecies(speciesData.species);
      } catch (err) {
        console.error("Failed to load species data:", err);
        setError("Failed to load species data. Please try again later.");
      } finally {
        setLoading(false);
      }
    };

    loadSpecies();
  }, []);

  // Set default species once the species list is loaded
  useEffect(() => {
    if (species.length > 0 && !formData.species) {
      setFormData((prev) => ({ ...prev, species: species[0] }));
    }
  }, [species, formData.species]);

  // Use the centralized validation hook
  const { errors, isValid } = useValidateForm(formData);

  const handleSubmit = (e) => {
    e.preventDefault();
    if (isValid && !isSubmitting) {
      onSubmit(formData);
    }
  };

  const updateFormData = (field, value) => {
    setFormData((prev) => ({ ...prev, [field]: value }));
  };

  const updateSequence = (index, field, value) => {
    setFormData((prev) => {
      const updatedSequences = [...prev.sequences];
      updatedSequences[index] = {
        ...updatedSequences[index],
        [field]: value,
      };
      return { ...prev, sequences: updatedSequences };
    });
  };

  const addSequence = () => {
    if (formData.numSequences < 10) {
      setFormData((prev) => ({
        ...prev,
        numSequences: prev.numSequences + 1,
        sequences: [
          ...prev.sequences,
          { sequence: "", primerName: "", mtkPart: "" },
        ],
      }));
    }
  };

  const removeSequence = () => {
    if (formData.numSequences > 1) {
      setFormData((prev) => ({
        ...prev,
        numSequences: prev.numSequences - 1,
        sequences: prev.sequences.slice(0, -1),
      }));
    }
  };

  // Reset form to defaults (while preserving species)
  const resetForm = () => {
    const defaults = getDefaultFormData();
    setFormData({
      ...defaults,
      species: formData.species,
    });
  };

  return (
    <form id="primer-form" onSubmit={handleSubmit}>
      {loading && <div className="loading-overlay">Loading species data...</div>}
      {error && <div className="error-message">{error}</div>}

      {/* Settings Card */}
      <Settings
        show={showSettings}
        onClose={() => setShowSettings(false)}
        formData={formData}
        updateFormData={updateFormData}
        availableSpecies={species}
      />

      {/* Settings Toggle Button */}
      <button
        type="button"
        className="settings-toggle"
        onClick={() => setShowSettings(!showSettings)}
      >
        ⚙️ Settings
      </button>

      {/* Template Sequence Section */}
      <TemplateSequence
        value={formData.templateSequence}
        onChange={(value) => updateFormData("templateSequence", value)}
      />

      {/* Sequences to Domesticate Section */}
      <SequenceTabs
        sequences={formData.sequences}
        updateSequence={updateSequence}
        addSequence={addSequence}
        removeSequence={removeSequence}
      />

      {/* Optionally, render global error messages */}
      <div className="global-errors">
        {Object.entries(errors).map(([field, message]) => (
          <div key={field} className="error-message">
            {message}
          </div>
        ))}
      </div>

      {/* Form Buttons */}
      <div className="button-container">
        <button
          type="submit"
          className={`btn btn-primary ${isSubmitting ? "processing" : ""}`}
          disabled={!isValid || isSubmitting || loading}
        >
          {isSubmitting ? (
            <>
              <span className="spinner"></span>
              Processing...
            </>
          ) : (
            "Generate Protocol"
          )}
        </button>
        <button
          type="button"
          className="btn btn-warning"
          onClick={resetForm}
        >
          Clear Form
        </button>
      </div>
    </form>
  );
}

export default Form;
