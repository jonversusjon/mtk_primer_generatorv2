// hooks/useValidateForm.js
import { useState, useEffect } from "react";
import { validateDnaSequence } from "../utils/dnaUtils";

const useValidateForm = (formData) => {
  const [errors, setErrors] = useState({});
  const [advisories, setAdvisories] = useState({});

  useEffect(() => {
    const newErrors = {};
    const newAdvisories = {};

    // Validate species selection (required)
    if (!formData.species || formData.species.trim() === "") {
      newErrors.species = "Species is required.";
    }

    // Validate template sequence (optional field)
    if (formData.templateSequence && formData.templateSequence.trim() !== "") {
      const tempValidation = validateDnaSequence(
        formData.templateSequence,
        false,
        false
      );
      // If not valid and not advisory, it's an error.
      if (!tempValidation.isValid && !tempValidation.isAdvisory) {
        newErrors.templateSequence = tempValidation.message;
      } else if (tempValidation.isAdvisory) {
        newAdvisories.templateSequence = tempValidation.message;
      }
    }

    // Validate sequences array
    if (!formData.sequences || formData.sequences.length === 0) {
      newErrors.sequences = "At least one sequence is required.";
    } else {
      formData.sequences.forEach((seq, index) => {
        // Validate sequence: non-empty and passes custom DNA rules
        if (!seq.sequence || seq.sequence.trim() === "") {
          newErrors[`sequences[${index}].sequence`] =
            "Sequence cannot be empty.";
        } else {
          const seqValidation = validateDnaSequence(seq.sequence, true, true);
          if (!seqValidation.isValid && !seqValidation.isAdvisory) {
            newErrors[`sequences[${index}].sequence`] = seqValidation.message;
          } else if (seqValidation.isAdvisory) {
            newAdvisories[`sequences[${index}].sequence`] =
              seqValidation.message;
          }
        }
        // Validate primer name (required)
        if (!seq.primerName || seq.primerName.trim() === "") {
          newErrors[`sequences[${index}].primerName`] =
            "Primer name is required.";
        }
        // Validate MTK Part (required)
        if (!seq.mtkPart || seq.mtkPart.trim() === "") {
          newErrors[`sequences[${index}].mtkPart`] = "MTK Part is required.";
        }
      });
    }

    setErrors(newErrors);
    setAdvisories(newAdvisories);
  }, [formData]);

  // Overall form is valid if there are no required errors.
  const isValid = Object.keys(errors).length === 0;
  return { errors, advisories, isValid };
};

export default useValidateForm;
