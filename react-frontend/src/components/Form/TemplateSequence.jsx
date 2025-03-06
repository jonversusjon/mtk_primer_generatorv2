import React, { useState, useEffect } from "react";
import { validateDnaSequence } from "../../utils/dnaUtils";
import SequenceInput from "../Form/SequenceInput";
import "../../styles/TemplateSequence.css";

function TemplateSequence({ value, onChange }) {
  const [validationMessage, setValidationMessage] = useState("");
  const [charCount, setCharCount] = useState(0);

  useEffect(() => {
    // Update character count when value changes
    setCharCount(value.length);

    // Validate template sequence (optional field)
    if (value) {
      const validationResult = validateDnaSequence(value, false, false);
      if (!validationResult.isValid) {
        setValidationMessage(validationResult.message);
      } else {
        setValidationMessage("");
      }
    } else {
      setValidationMessage("");
    }
  }, [value]);

  const handleChange = (e) => {
    const newValue = e.target.value;
    onChange(newValue);
  };

  return (
    <div className="template-sequence-section">
      <div className="form-row sequence-header">
        <div className="section-title">
          <h3>Template Sequence</h3>
          <span className="optional-badge">(optional)</span>
        </div>

        <div className="char-count">Length: {charCount} bp</div>
      </div>

      <SequenceInput
        id="templateSequence"
        value={value}
        onChange={handleChange}
        placeholder="Paste your template DNA sequence here (optional)"
      />

      {validationMessage && (
        <div className="validation-message error">{validationMessage}</div>
      )}
    </div>
  );
}

export default TemplateSequence;
