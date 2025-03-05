import React, { useState, useEffect } from "react";
import SequenceInput from "../Form/SequenceInput";
import { validateDnaSequence } from "../../utils/dnaUtils";
import "../../styles/SequenceTab.css";

function SequenceTab({ sequence, index, updateSequence, mtkPartOptions }) {
  // Instead of just a string, we store both the message and whether it's advisory
  const [validation, setValidation] = useState({
    message: "",
    isAdvisory: false,
  });
  const [charCount, setCharCount] = useState(0);
  const [useSeparateParts, setUseSeparateParts] = useState(false);

  // Initialize values on component mount
  useEffect(() => {
    // If mtkPartLeft or mtkPartRight doesn't exist (first load), initialize them
    if (!sequence.mtkPartLeft || !sequence.mtkPartRight) {
      // Initialize with empty strings or use a reasonable default from options
      const initialValue = mtkPartOptions.length > 0 ? mtkPartOptions[0] : "";
      updateSequence("mtkPartLeft", initialValue);
      updateSequence("mtkPartRight", initialValue);
    }
  }, [
    sequence.mtkPartLeft,
    sequence.mtkPartRight,
    updateSequence,
    index,
    mtkPartOptions,
  ]);

  useEffect(() => {
    // Update character count
    setCharCount(sequence.sequence.length);

    // Validate DNA sequence if one is provided
    if (sequence.sequence) {
      const validationResult = validateDnaSequence(sequence.sequence, true); // true for required validation
      if (validationResult.isValid) {
        setValidation({ message: "", isAdvisory: false });
      } else {
        // Store the message and its advisory status
        setValidation({
          message: validationResult.message,
          isAdvisory: validationResult.isAdvisory,
        });
      }
    } else {
      setValidation({ message: "", isAdvisory: false });
    }
  }, [sequence.sequence]);

  const handleSequenceChange = (e) => {
    updateSequence("sequence", e.target.value);
  };

  const handlePrimerNameChange = (e) => {
    updateSequence("primerName", e.target.value);
  };

  const handleMtkPartChange = (partType) => (e) => {
    const newValue = e.target.value;

    if (useSeparateParts) {
      // Only update the specific part that changed
      updateSequence(partType, newValue);
    } else {
      // In single mode, update both left and right to maintain synchronization
      updateSequence("mtkPartLeft", newValue);
      updateSequence("mtkPartRight", newValue);
    }
  };

  const handleToggleChange = (e) => {
    setUseSeparateParts(e.target.checked);
  };

  return (
    <div className="sequence-tab-content">
      <div className="sequence-header">
        <label htmlFor={`sequence-${index}`}>Sequence {index + 1}:</label>
        <div
          className="char-count"
          style={{ display: charCount > 0 ? "inline" : "none" }}
        >
          Length: {charCount} bp
        </div>
      </div>

      <SequenceInput
        id={`sequence-${index}`}
        value={sequence.sequence}
        onChange={handleSequenceChange}
        placeholder="Paste your DNA sequence here"
      />

      {/* Display the message with a conditional className */}
      {validation.message && (
        <div
          className={`validation-message ${
            validation.isAdvisory ? "advisory" : "error"
          }`}
        >
          {validation.message}
        </div>
      )}

      <div className="sequence-metadata">
        <div className="form-group">
          <label htmlFor={`primer-name-${index}`}>Primer Name:</label>
          <input
            type="text"
            id={`primer-name-${index}`}
            value={sequence.primerName}
            onChange={handlePrimerNameChange}
            className="form-control"
            placeholder="Enter primer name"
          />
        </div>

        <div className="form-group toggle-container">
          <label className="toggle-label">
            <input
              type="checkbox"
              checked={useSeparateParts}
              onChange={handleToggleChange}
            />
            Use separate left/right part numbers
          </label>
        </div>

        {!useSeparateParts ? (
          // Single MTK part number selector (but updates both left/right internally)
          <div className="form-group">
            <label htmlFor={`mtk-part-${index}`}>MTK Part Number:</label>
            <select
              id={`mtk-part-${index}`}
              value={sequence.mtkPartLeft || ""}
              onChange={handleMtkPartChange("mtkPartLeft")}
              className="form-control"
            >
              {mtkPartOptions.map((part) => (
                <option key={part} value={part}>
                  {part}
                </option>
              ))}
            </select>
          </div>
        ) : (
          // Two MTK part number selectors
          <>
            <div className="form-group">
              <label htmlFor={`mtk-part-left-${index}`}>
                MTK Part Number Left:
              </label>
              <select
                id={`mtk-part-left-${index}`}
                value={sequence.mtkPartLeft || ""}
                onChange={handleMtkPartChange("mtkPartLeft")}
                className="form-control"
              >
                {mtkPartOptions.map((part) => (
                  <option key={part} value={part}>
                    {part}
                  </option>
                ))}
              </select>
            </div>
            <div className="form-group">
              <label htmlFor={`mtk-part-right-${index}`}>
                MTK Part Number Right:
              </label>
              <select
                id={`mtk-part-right-${index}`}
                value={sequence.mtkPartRight || ""}
                onChange={handleMtkPartChange("mtkPartRight")}
                className="form-control"
              >
                {mtkPartOptions.map((part) => (
                  <option key={part} value={part}>
                    {part}
                  </option>
                ))}
              </select>
            </div>
          </>
        )}
      </div>
    </div>
  );
}

export default SequenceTab;
