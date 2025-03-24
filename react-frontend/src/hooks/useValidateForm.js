import { useState, useEffect } from "react";
import { validateDnaSequence } from "../utils/dnaUtils";

const useValidateForm = (formData, shouldValidate = true) => {
  const [errors, setErrors] = useState({});
  const [advisories, setAdvisories] = useState({});

  useEffect(() => {

    if (!shouldValidate) {
      console.log("Validation skipped because shouldValidate is false.");
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
        console.log(
          "Validation error: species is empty. formData.species:",
          formData.species
        );
      }

      // Validate template sequence (optional field)
      if (
        formData.templateSequence &&
        formData.templateSequence.trim() !== ""
      ) {
        const tempValidation = validateDnaSequence(
          formData.templateSequence,
          false,
          false
        );
        if (!tempValidation.isValid && !tempValidation.isAdvisory) {
          newErrors.templateSequence = tempValidation.message;
          console.log(
            "Validation error in templateSequence:",
            formData.templateSequence,
            tempValidation.message
          );
        } else if (tempValidation.isAdvisory) {
          newAdvisories.templateSequence = tempValidation.message;
          console.log(
            "Advisory for templateSequence:",
            formData.templateSequence,
            tempValidation.message
          );
        }
      } else {
        console.log(
          "Template sequence is empty or not provided; skipping validation."
        );
      }

      // Validate sequencesToDomesticate array
      if (
        !formData.sequencesToDomesticate ||
        formData.sequencesToDomesticate.length === 0
      ) {
        newErrors.sequencesToDomesticate = "At least one sequence is required.";
        console.log("Validation error: sequencesToDomesticate array is empty.");
      } else {
        formData.sequencesToDomesticate.forEach((seq, index) => {
          const trimmedSequence =
            typeof seq.sequence === "string" ? seq.sequence.trim() : "";
          if (!trimmedSequence) {
            newErrors[`sequencesToDomesticate[${index}].sequence`] =
              "Sequence cannot be empty.";
            console.log(
              `Validation error at sequence index ${index}: sequence is empty.`
            );
          } else {
            const seqValidation = validateDnaSequence(seq.sequence, true, true);

            if (!seqValidation.isValid && !seqValidation.isAdvisory) {
              newErrors[`sequencesToDomesticate[${index}].sequence`] =
                seqValidation.message;
              console.log(
                `Validation error at sequence index ${index}:`,
                seqValidation.message
              );
            } else if (seqValidation.isAdvisory) {
              newAdvisories[`sequencesToDomesticate[${index}].sequence`] =
                seqValidation.message;
              console.log(
                `Validation advisory at sequence index ${index}:`,
                seqValidation.message
              );
            } 
          }

          // Validate primer name (required)
          if (!seq.primerName || seq.primerName.trim() === "") {
            newErrors[`sequencesToDomesticate[${index}].primerName`] =
              "Primer name is required.";
            console.log(
              `Validation error at sequence index ${index}: Primer name is empty.`
            );
          } 

          // Validate MTK Part Left (required)
          if (!seq.mtkPartLeft || seq.mtkPartLeft.trim() === "") {
            newErrors[`sequencesToDomesticate[${index}].mtkPartLeft`] =
              "MTK Part Left is required.";
            console.log(
              `Validation error at sequence index ${index}: MTK Part Left is empty.`
            );
          }

          // Validate MTK Part Right (required)
          if (!seq.mtkPartRight || seq.mtkPartRight.trim() === "") {
            newErrors[`sequencesToDomesticate[${index}].mtkPartRight`] =
              "MTK Part Right is required.";
            console.log(
              `Validation error at sequence index ${index}: MTK Part Right is empty.`
            );
          }
        });
      }

      console.log("Final validation errors:", newErrors);
      setErrors(newErrors);
      setAdvisories(newAdvisories);
    }, 300);

    return () => clearTimeout(timer);
  }, [formData, shouldValidate]);

  // Overall form is valid if there are no required errors.
  const isValid = Object.keys(errors).length === 0;
  if (!isValid) {
    console.log("Form validation errors:", errors);
  }
  return { errors, advisories, isValid };
};

export default useValidateForm;
