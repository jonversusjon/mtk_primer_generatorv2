import React, { useState, useEffect, useRef, useCallback } from "react";
import TemplateSequence from "./TemplateSequence";
import SequenceTabs from "./SequenceTabs";
import { fetchAvailableSpecies } from "../../api/api";

const Form = ({
  onSubmit,
  formData,
  setFormData,
  errors,
  isValid,
  activeTabIndex,
  setActiveTabIndex,
}) => {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [species, setSpecies] = useState([]);
  const speciesDefaultSet = useRef(false);

  // Load species and update formData with availableSpecies
  useEffect(() => {
    const loadSpecies = async () => {
      setLoading(true);
      setError(null);
      try {
        const speciesData = await fetchAvailableSpecies();
        setSpecies(speciesData.species);
        // Update formData so that Settings (in Sidebar) can use the species list
        setFormData((prev) => ({
          ...prev,
          availableSpecies: speciesData.species,
        }));
      } catch (err) {
        console.error("Failed to load species data:", err);
        setError("Failed to load species data. Please try again later.");
      } finally {
        setLoading(false);
      }
    };
    loadSpecies();
  }, [setFormData]);

  // Set default species if not already set
  useEffect(() => {
    if (
      species.length > 0 &&
      (!formData.species || formData.species === "") &&
      !speciesDefaultSet.current
    ) {
      setFormData((prev) => ({ ...prev, species: species[0] }));
      speciesDefaultSet.current = true;
    }
  }, [species, formData.species, setFormData]);

  const handleSubmit = (e) => {
    e.preventDefault();
    // Ensure species is set in the submitted data
    const finalFormData =
      (!formData.species || formData.species === "") && species.length > 0
        ? { ...formData, species: species[0] }
        : formData;
    if (isValid) {
      onSubmit(finalFormData);
    }
  };

  // Functional updates to formData
  const updateFormData = useCallback(
    (field, value) => {
      setFormData((prev) => ({ ...prev, [field]: value }));
    },
    [setFormData]
  );

  const updateSequence = useCallback(
    (index, field, value) => {
      setFormData((prev) => {
        const updatedSequences = [...prev.sequencesToDomesticate];
        if (updatedSequences[index]?.[field] === value) return prev;
        updatedSequences[index] = {
          ...updatedSequences[index],
          [field]: value,
        };
        return { ...prev, sequencesToDomesticate: updatedSequences };
      });
    },
    [setFormData]
  );

  const addSequence = useCallback(() => {
    setFormData((prev) => ({
      ...prev,
      sequencesToDomesticate: [
        ...prev.sequencesToDomesticate,
        { sequence: "", primerName: "", mtkPartLeft: "", mtkPartRight: "" },
      ],
    }));
  }, [setFormData]);

  const removeSequence = useCallback(() => {
    setFormData((prev) => {
      if (prev.sequencesToDomesticate.length <= 1) return prev;
      return {
        ...prev,
        sequencesToDomesticate: prev.sequencesToDomesticate.slice(0, -1),
      };
    });
  }, [setFormData]);

  const resetForm = useCallback(() => {
    setFormData((prev) => ({ ...prev, templateSequence: "", species: "" }));
  }, [setFormData]);

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
};

export default Form;
