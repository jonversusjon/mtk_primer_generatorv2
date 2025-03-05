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

    // Validate sequencesToDomesticate array
    if (
      !formData.sequencesToDomesticate ||
      formData.sequencesToDomesticate.length === 0
    ) {
      newErrors.sequencesToDomesticate = "At least one sequence is required.";
    } else {
      formData.sequencesToDomesticate.forEach((seq, index) => {
        // Validate sequence: non-empty and passes custom DNA rules
        if (!seq.sequence || seq.sequence.trim() === "") {
          newErrors[`sequencesToDomesticate[${index}].sequence`] =
            "Sequence cannot be empty.";
        } else {
          const seqValidation = validateDnaSequence(seq.sequence, true, true);
          if (!seqValidation.isValid && !seqValidation.isAdvisory) {
            newErrors[`sequencesToDomesticate[${index}].sequence`] =
              seqValidation.message;
          } else if (seqValidation.isAdvisory) {
            newAdvisories[`sequencesToDomesticate[${index}].sequence`] =
              seqValidation.message;
          }
        }

        // Validate primer name (required)
        if (!seq.primerName || seq.primerName.trim() === "") {
          newErrors[`sequencesToDomesticate[${index}].primerName`] =
            "Primer name is required.";
        }

        // Validate MTK Part Left (required)
        if (!seq.mtkPartLeft || seq.mtkPartLeft.trim() === "") {
          newErrors[`sequencesToDomesticate[${index}].mtkPartLeft`] =
            "MTK Part Left is required.";
        }

        // Validate MTK Part Right (required)
        if (!seq.mtkPartRight || seq.mtkPartRight.trim() === "") {
          newErrors[`sequencesToDomesticate[${index}].mtkPartRight`] =
            "MTK Part Right is required.";
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
