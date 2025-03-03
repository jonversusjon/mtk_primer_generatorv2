import React, { useState, useEffect } from "react";
import SequenceInput from "../Form/SequenceInput";
import { validateDnaSequence } from "../../utils/dnaUtils";
import "../../styles/SequenceTab.css";

function SequenceTab({ sequence, index, updateSequence, mtkPartOptions }) {
  const [validationMessage, setValidationMessage] = useState("");
  const [charCount, setCharCount] = useState(0);

  useEffect(() => {
    // Update character count
    setCharCount(sequence.sequence.length);

    // Validate DNA sequence
    if (sequence.sequence) {
      const validationResult = validateDnaSequence(sequence.sequence, true); // true for mandatory validation
      setValidationMessage(
        validationResult.isValid ? "" : validationResult.message
      );
    } else {
      setValidationMessage("");
    }
  }, [sequence.sequence]);

  const handleSequenceChange = (e) => {
    updateSequence("sequence", e.target.value);
  };

  const handlePrimerNameChange = (e) => {
    updateSequence("primerName", e.target.value);
  };

  const handleMtkPartChange = (e) => {
    updateSequence("mtkPart", e.target.value);
  };

  return (
    <>
      <div className="sequence-header">
        <label htmlFor={`sequence${index}`} className="form-label">
          Sequence {index + 1}:
        </label>
        <span
          id={`charCount${index}`}
          className="char-count-label"
          style={{ display: charCount > 0 ? "inline" : "none" }}
        >
          Length: {charCount} bp
        </span>
      </div>

      <SequenceInput
        id={`sequence${index}`}
        name={`sequences[${index}][sequence]`}
        value={sequence.sequence}
        onChange={handleSequenceChange}
        placeholder="Paste DNA sequence..."
        className="dynamic-sequence-input"
        required
      />

      {validationMessage && (
        <div className="validation-feedback invalid">{validationMessage}</div>
      )}

      <div className="form-group">
        <label htmlFor={`primerName${index}`}>Primer Name:</label>
        <input
          type="text"
          id={`primerName${index}`}
          name={`sequences[${index}][primerName]`}
          placeholder="Default Primer Name"
          value={sequence.primerName}
          onChange={handlePrimerNameChange}
          className="primer-name-input"
          required
        />

        <label htmlFor={`mtkPart${index}`}>MTK Part Number:</label>
        <select
          id={`mtkPart${index}`}
          name={`sequences[${index}][mtkPart]`}
          value={sequence.mtkPart}
          onChange={handleMtkPartChange}
          className="mtk-part-select"
          required
        >
          {mtkPartOptions.map((part) => (
            <option key={part} value={part}>
              {part}
            </option>
          ))}
        </select>
      </div>
    </>
  );
}

export default SequenceTab;
