import React, { useCallback } from "react";
import TemplateSequence from "./TemplateSequence";
import SequenceTabs from "./SequenceTabs";

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
  // No local species state or fetching here.

  const handleSubmit = (e) => {
    e.preventDefault();
    // Ensure species is set before submitting.
    // If formData.species is empty, default to the first species in availableSpecies.
    const finalFormData =
      (!formData.species || formData.species.trim() === "") &&
      formData.availableSpecies.length > 0
        ? { ...formData, species: formData.availableSpecies[0] }
        : formData;
    if (isValid && initialized) {
      onSubmit(finalFormData);
    }
  };

  // Curried updateSequence function so that SequenceTabs (and its children) can call it as updateSequence(index)(field, value)
  const updateSequence = useCallback(
    (index) => (field, value) => {
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
      {!initialized && (
        <div className="loading-overlay">Loading required data...</div>
      )}

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
          disabled={!isValid || !initialized}
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
