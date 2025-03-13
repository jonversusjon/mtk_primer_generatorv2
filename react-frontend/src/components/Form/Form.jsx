import React, { useState, useEffect, useRef, useCallback } from "react";
import TemplateSequence from "./TemplateSequence";
import SequenceTabs from "./SequenceTabs";
import { fetchAvailableSpecies } from "../../api/api";

function Form({
  onSubmit,
  formData,
  setFormData,
  errors,
  isValid,
  activeTabIndex,
  setActiveTabIndex,
}) {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [species, setSpecies] = useState([]);
  const speciesDefaultSet = useRef(false);

  // Load species
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

  // Ensure default species is set once
  useEffect(() => {
    if (
      species.length > 0 &&
      (!formData.species || formData.species === "") &&
      !speciesDefaultSet.current
    ) {
      setFormData({ ...formData, species: species[0] });
      speciesDefaultSet.current = true;
    }
  }, [species, formData, setFormData]);

  const handleSubmit = (e) => {
    e.preventDefault();
    // Safeguard for species
    const finalFormData =
      (!formData.species || formData.species === "") && species.length > 0
        ? { ...formData, species: species[0] }
        : formData;

    if (isValid) {
      onSubmit(finalFormData);
    }
  };

  // Wrap updates in useCallback if desired
  const updateFormData = useCallback(
    (field, value) => {
      setFormData({ ...formData, [field]: value });
    },
    [formData, setFormData]
  );

  const updateSequence = useCallback(
    (index, field, value) => {
      if (formData.sequencesToDomesticate[index]?.[field] === value) return;
      const updatedSequences = [...formData.sequencesToDomesticate];
      updatedSequences[index] = {
        ...updatedSequences[index],
        [field]: value,
      };
      setFormData({ ...formData, sequencesToDomesticate: updatedSequences });
    },
    [formData, setFormData]
  );

  const addSequence = useCallback(() => {
    const updatedSequences = [
      ...formData.sequencesToDomesticate,
      { sequence: "", primerName: "", mtkPartLeft: "", mtkPartRight: "" },
    ];
    setFormData({ ...formData, sequencesToDomesticate: updatedSequences });
  }, [formData, setFormData]);

  const removeSequence = useCallback(() => {
    if (formData.sequencesToDomesticate.length > 1) {
      setFormData({
        ...formData,
        sequencesToDomesticate: formData.sequencesToDomesticate.slice(0, -1),
      });
    }
  }, [formData, setFormData]);

  const resetForm = useCallback(() => {
    // If you have "default form data" you can pass that in
    setFormData({ ...formData, templateSequence: "", species: "" });
  }, [formData, setFormData]);

  return (
    <form id="primer-form" onSubmit={handleSubmit}>
      {loading && (
        <div className="loading-overlay">Loading species data...</div>
      )}
      {error && <div className="error-message">{error}</div>}

      <TemplateSequence
        value={formData.templateSequence}
        onChange={(value) => updateFormData("templateSequence", value)}
      />

      <SequenceTabs
        sequencesToDomesticate={formData.sequencesToDomesticate}
        updateSequence={updateSequence}
        addSequence={addSequence}
        removeSequence={removeSequence}
        activeTabIndex={activeTabIndex}
        onTabChange={setActiveTabIndex}
      />

      <div className="button-container">
        <button
          type="submit"
          className="btn btn-primary"
          disabled={!isValid || loading}
        >
          Generate Protocol
        </button>
        <button type="button" className="btn btn-warning" onClick={resetForm}>
          Clear Form
        </button>
      </div>
    </form>
  );
}

export default Form;
