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
      const validationResult = validateDnaSequence(value);
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
    <div className="form-sub-container">
      <div className="section-header">
        <div className="template-sequence-header">
          <h2 className="form-sub-container-title">Template Sequence</h2>
          <h4>(optional)</h4>
        </div>
      </div>

      <div className="paste-template-sequence-header">
        <label htmlFor="templateSequence" className="form-label">
          Paste Template Sequence:
        </label>
        <span id="charCount" className="char-count-label">
          Length: {charCount} bp
        </span>
      </div>

      <SequenceInput
        id="templateSequence"
        name="templateSequence"
        value={value}
        onChange={handleChange}
        placeholder="Enter DNA sequence here..."
        className="sequence-input"
      />

      {validationMessage && (
        <div className="error-message" data-field="templateSequence">
          {validationMessage}
        </div>
      )}
    </div>
  );
}

export default TemplateSequence;
