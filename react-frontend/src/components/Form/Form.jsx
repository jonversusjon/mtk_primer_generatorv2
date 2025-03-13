// Form.jsx
import React, { useState, useEffect, useRef } from "react";
import TemplateSequence from "./TemplateSequence";
import SequenceTabs from "./SequenceTabs";
import Settings from "./Settings";
import { fetchAvailableSpecies } from "../../api/api";
import useValidateForm from "../../hooks/useValidateForm";

// formSchema.js - Define your data structure clearly in one place
export const formSchema = {
  templateSequence: {
    type: "string",
    description: "The template DNA sequence",
  },
  species: { type: "string", description: "Species selection" },
  kozak: { type: "string", description: "Kozak sequence type" },
  max_mut_per_site: {
    type: "number",
    description: "Maximum mutations per site",
  },
  verbose_mode: { type: "boolean", description: "Enable verbose output" },
  sequencesToDomesticate: {
    type: "array",
    items: {
      sequence: { type: "string", description: "DNA sequence" },
      primerName: { type: "string", description: "Name of the primer" },
      mtkPartLeft: { type: "string", description: "Left MTK part number" },
      mtkPartRight: { type: "string", description: "Right MTK part number" },
    },
  },
};

function Form({
  onSubmit,
  isSubmitting,
  formData,
  onChange,
  showSettings,
  setShowSettings,
}) {
  // State for species loading and errors
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [species, setSpecies] = useState([]);

  // Ref to ensure we set default species only once
  const speciesDefaultSet = useRef(false);

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

  // Set default species in formData only once if it's blank
  useEffect(() => {
    if (
      species.length > 0 &&
      (!formData.species || formData.species === "") &&
      !speciesDefaultSet.current
    ) {
      onChange({ ...formData, species: species[0] });
      speciesDefaultSet.current = true;
    }
  }, [species, formData, onChange]);

  // Use the centralized validation hook
  const { errors, isValid } = useValidateForm(formData);

  const handleSubmit = (e) => {
    e.preventDefault();

    // Ensure species is set before submission.
    // This is an extra safeguard so the backend never gets a blank species.
    const finalFormData =
      (!formData.species || formData.species === "") && species.length > 0
        ? { ...formData, species: species[0] }
        : formData;

    if (isValid && !isSubmitting) {
      console.log("Form data being sent:", finalFormData);
      onSubmit(finalFormData);
    }
  };

  const updateFormData = (field, value) => {
    onChange({ ...formData, [field]: value });
  };

  const updateSequence = (index, field, value) => {
    // Prevent unnecessary updates if the value hasn't changed
    if (formData.sequencesToDomesticate[index]?.[field] === value) {
      return;
    }
    const updatedSequences = [...formData.sequencesToDomesticate];
    updatedSequences[index] = { ...updatedSequences[index], [field]: value };
    onChange({ ...formData, sequencesToDomesticate: updatedSequences });
  };

  const addSequence = () => {
    const updatedSequences = [
      ...formData.sequencesToDomesticate,
      { sequence: "", primerName: "", mtkPartLeft: "", mtkPartRight: "" },
    ];
    onChange({ ...formData, sequencesToDomesticate: updatedSequences });
  };

  const removeSequence = () => {
    if (formData.sequencesToDomesticate.length > 1) {
      onChange({
        ...formData,
        sequencesToDomesticate: formData.sequencesToDomesticate.slice(0, -1),
      });
    }
  };

  // Reset form to defaults (while preserving species)
  const resetForm = () => {
    onChange({ ...formData });
  };

  return (
    <form id="primer-form" onSubmit={handleSubmit}>
      <div className="form-header">
        <h2 className="primer-form-title">Primer Design Form</h2>
        <button
          type="button"
          className="settings-toggle"
          onClick={() => setShowSettings(!showSettings)}
          aria-label={showSettings ? "Close settings" : "Open settings"}
        >
          <span className="icon">⚙️</span>
          <span className="button-text">Settings</span>
        </button>
      </div>

      {loading && (
        <div className="loading-overlay">Loading species data...</div>
      )}
      {error && <div className="error-message">{error}</div>}

      {/* Settings Card */}
      <Settings
        show={showSettings}
        onClose={() => setShowSettings(false)}
        formData={formData}
        updateFormData={updateFormData}
        availableSpecies={species}
      />

      {/* Template Sequence Section */}
      <TemplateSequence
        value={formData.templateSequence}
        onChange={(value) => updateFormData("templateSequence", value)}
      />

      {/* Sequences to Domesticate Section */}
      <SequenceTabs
        sequencesToDomesticate={formData.sequencesToDomesticate}
        updateSequence={updateSequence}
        addSequence={addSequence}
        removeSequence={removeSequence}
      />

      {/* Optionally, render global error messages */}
      {!loading && species.length > 0 && (
        <div className="global-errors">
          {Object.entries(errors).map(([field, message]) => (
            <div key={field} className="error-message">
              {message}
            </div>
          ))}
        </div>
      )}

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
        <button type="button" className="btn btn-warning" onClick={resetForm}>
          Clear Form
        </button>
      </div>
    </form>
  );
}

export default Form;
