import { useState, useEffect } from "react";
import { validateDnaSequence } from "../utils/dnaUtils";

const useValidateForm = (formData, shouldValidate = true) => {
  const [errors, setErrors] = useState({});
  const [advisories, setAdvisories] = useState({});

  useEffect(() => {
    console.log("useValidateForm - formData:", formData, "shouldValidate:", shouldValidate);

    if (!shouldValidate) {
      // Clear errors while form is still loading or initializing
      setErrors({});
      setAdvisories({});
      return;
    }

    // Delay validation slightly to accommodate asynchronous initialization
    const timer = setTimeout(() => {
      const newErrors = {};
      const newAdvisories = {};

      // Validate species selection (required)
      if (!formData.species || formData.species.trim() === "") {
        newErrors.species = "Species is required.";
      }

      // Validate template sequence (optional field)
      if (formData.templateSequence && formData.templateSequence.trim() !== "") {
        const tempValidation = validateDnaSequence(formData.templateSequence, false, false);
        if (!tempValidation.isValid && !tempValidation.isAdvisory) {
          newErrors.templateSequence = tempValidation.message;
        } else if (tempValidation.isAdvisory) {
          newAdvisories.templateSequence = tempValidation.message;
        }
      }

      // Validate sequencesToDomesticate array
      if (!formData.sequencesToDomesticate || formData.sequencesToDomesticate.length === 0) {
        newErrors.sequencesToDomesticate = "At least one sequence is required.";
      } else {
        formData.sequencesToDomesticate.forEach((seq, index) => {
          const trimmedSequence = typeof seq.sequence === "string" ? seq.sequence.trim() : "";
          if (!trimmedSequence) {
            newErrors[`sequencesToDomesticate[${index}].sequence`] = "Sequence cannot be empty.";
          } else {
            const seqValidation = validateDnaSequence(seq.sequence, true, true);
            console.log(`Validation for sequence at index ${index}:`, seqValidation);
            if (!seqValidation.isValid && !seqValidation.isAdvisory) {
              newErrors[`sequencesToDomesticate[${index}].sequence`] = seqValidation.message;
            } else if (seqValidation.isAdvisory) {
              newAdvisories[`sequencesToDomesticate[${index}].sequence`] = seqValidation.message;
            }
          }

          // Validate primer name (required)
          if (!seq.primerName || seq.primerName.trim() === "") {
            newErrors[`sequencesToDomesticate[${index}].primerName`] = "Primer name is required.";
          }

          // Validate MTK Part Left (required)
          if (!seq.mtkPartLeft || seq.mtkPartLeft.trim() === "") {
            newErrors[`sequencesToDomesticate[${index}].mtkPartLeft`] = "MTK Part Left is required.";
          }

          // Validate MTK Part Right (required)
          if (!seq.mtkPartRight || seq.mtkPartRight.trim() === "") {
            newErrors[`sequencesToDomesticate[${index}].mtkPartRight`] = "MTK Part Right is required.";
          }
        });
      }

      console.log("Validation errors:", newErrors);
      setErrors(newErrors);
      setAdvisories(newAdvisories);
    }, 300);

    return () => clearTimeout(timer);
  }, [formData, shouldValidate]);

  // Overall form is valid if there are no required errors.
  const isValid = Object.keys(errors).length === 0;
  return { errors, advisories, isValid };
};

export default useValidateForm;
