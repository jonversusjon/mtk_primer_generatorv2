import React, { useState, useEffect, useCallback } from "react";
import TemplateSequence from "./TemplateSequence";
import SequenceTabs from "./SequenceTabs";
import { fetchAvailableSpecies } from "../../api/api";

const Form = ({
  onSubmit,
  formData,
  setFormData,
  errors,
  isValid,
  initialized,
  activeTabIndex,
  setActiveTabIndex,
}) => {
  const [loadingSpecies, setLoadingSpecies] = useState(false);
  const [speciesError, setSpeciesError] = useState(null);
  const [species, setSpecies] = useState([]);

  // Load species and update formData with availableSpecies
  useEffect(() => {
    const loadSpecies = async () => {
      setLoadingSpecies(true);
      setSpeciesError(null);
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
        setSpeciesError("Failed to load species data. Please try again later.");
      } finally {
        setLoadingSpecies(false);
      }
    };
    loadSpecies();
  }, [setFormData]);

  // Set default species if available and not already set
  useEffect(() => {
    if (species.length > 0 && (!formData.species || formData.species.trim() === "")) {
      setFormData((prev) => ({ ...prev, species: species[0] }));
    }
  }, [species, formData.species, setFormData]);

  const handleSubmit = (e) => {
    e.preventDefault();
    // Ensure species is set before submitting
    const finalFormData =
      (!formData.species || formData.species.trim() === "") && species.length > 0
        ? { ...formData, species: species[0] }
        : formData;
    // Only submit if the form is valid and fully initialized
    if (isValid && initialized) {
      onSubmit(finalFormData);
    }
  };

  // Curried updateSequence function so that SequenceTabs (and its children) can call it as updateSequence(index)(field, value)
  const updateSequence = useCallback(
    (index) => (field, value) => {
      console.log(`updateSequence - index: ${index}, field: ${field}, value: ${value}`);
      setFormData((prev) => {
        const updatedSequences = [...prev.sequencesToDomesticate];
        if (updatedSequences[index]?.[field] === value) return prev;
        updatedSequences[index] = {
          ...updatedSequences[index],
          [field]: value,
        };
        console.log("Updated formData:", { ...prev, sequencesToDomesticate: updatedSequences });
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
      {(loadingSpecies || !initialized) && (
        <div className="loading-overlay">Loading required data...</div>
      )}
      {speciesError && <div className="error-message">{speciesError}</div>}

      <TemplateSequence
        value={formData.templateSequence}
        onChange={(value) =>
          setFormData((prev) => ({ ...prev, templateSequence: value }))
        }
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
          disabled={!isValid || loadingSpecies || !initialized}
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
