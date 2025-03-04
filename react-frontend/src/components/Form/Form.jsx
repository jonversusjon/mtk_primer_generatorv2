import React, { useState, useEffect, useCallback } from "react";
import TemplateSequence from "./TemplateSequence";
import SequenceTabs from "./SequenceTabs";
import Settings from "./Settings";
import { fetchAvailableSpecies } from "../../api/api";
import { defaultParameters } from "../../config/defaultParameters";
import "../../styles/form.css";

function Form({ onSubmit, isSubmitting }) {
  // State variables for species loading
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [species, setSpecies] = useState([]);

  // Custom hook for loading defaults based on environment mode
  const useLoadDefaults = () => {
    return useCallback(() => {
      console.log("Loading defaults with parameters:", defaultParameters);

      // Map the default sequences from config
      const initialSequences = [];

      // If we have test sequences (in TESTING mode), populate them
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
        species: "", // Will be set after fetching available species
        kozak: "MTK",
        max_mut_per_site: 1,
        verbose_mode: true,
        sequences: initialSequences,
      };
    }, []);
  };

  const loadDefaults = useLoadDefaults();
  console.log("Species:", species);

  // Form state initialization with default parameters
  const [formData, setFormData] = useState(loadDefaults());

  // Settings state
  const [showSettings, setShowSettings] = useState(false);
  const [availableSpecies, setAvailableSpecies] = useState([]);

  // Form validation
  const [isFormValid, setIsFormValid] = useState(false);

  // Load available species when component mounts
  useEffect(() => {
    const loadSpecies = async () => {
      try {
        setLoading(true);
        setError(null);
        const speciesData = await fetchAvailableSpecies();
        setSpecies(speciesData);
        setAvailableSpecies(speciesData); // Update availableSpecies for Settings component
      } catch (error) {
        console.error("Failed to load species data:", error);
        setError("Failed to load species data. Please try again later.");
      } finally {
        setLoading(false);
      }
    };

    loadSpecies();
  }, []);

  // Memoize validateForm function to prevent unnecessary rerenders
  const validateForm = useCallback(() => {
    // Check if template sequence is valid (optional)
    // Check if at least one sequence is provided and valid
    // Check if species is selected
    const hasValidSequence = formData.sequences.some(
      (seq) => seq.sequence && seq.primerName && seq.mtkPart
    );

    const isValid = formData.species && hasValidSequence;
    setIsFormValid(isValid);
  }, [formData]);

  // Update form validity when formData changes
  useEffect(() => {
    validateForm();
  }, [validateForm]);

  const handleSubmit = (e) => {
    e.preventDefault();
    if (isFormValid && !isSubmitting) {
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

  // Reset form to defaults
  const resetForm = () => {
    const defaults = loadDefaults();
    // Preserve the current species selection
    setFormData({
      ...defaults,
      species: formData.species,
    });
  };

  return (
    <form id="primer-form" onSubmit={handleSubmit}>
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
        availableSpecies={availableSpecies}
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

      {/* Form Buttons */}
      <div className="button-container">
        <button
          type="submit"
          className={`btn btn-primary ${isSubmitting ? "processing" : ""}`}
          id="runDesignPrimerBtn"
          disabled={!isFormValid || isSubmitting || loading}
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
          id="clearForm"
          onClick={resetForm}
        >
          Clear Form
        </button>
      </div>
    </form>
  );
}

export default Form;
