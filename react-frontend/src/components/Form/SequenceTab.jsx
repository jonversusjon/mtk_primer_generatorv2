import React, { useState, useEffect, useRef } from "react";
import SequenceInput from "../Form/SequenceInput";
import { validateDnaSequence } from "../../utils/dnaUtils";
import "../../styles/SequenceTab.css";

function SequenceTab({ sequence, index, updateSequence, mtkPartOptions }) {
  const [validation, setValidation] = useState({
    message: "",
    isAdvisory: false,
  });
  // Instead of a state for charCount, compute it directly:
  const charCount = sequence.sequence ? sequence.sequence.length : 0;
  const [useSeparateParts, setUseSeparateParts] = useState(false);

  const hasInitialized = useRef(false);
  useEffect(() => {
    if (!hasInitialized.current && mtkPartOptions.length > 0) {
      // Only update if the current values are empty
      if (!sequence.mtkPartLeft || sequence.mtkPartLeft.trim() === "") {
        updateSequence("mtkPartLeft", mtkPartOptions[0]);
      }
      if (!sequence.mtkPartRight || sequence.mtkPartRight.trim() === "") {
        updateSequence("mtkPartRight", mtkPartOptions[0]);
      }
      hasInitialized.current = true;
    }
  }, [
    mtkPartOptions,
    updateSequence,
    sequence.mtkPartLeft,
    sequence.mtkPartRight,
  ]);

  useEffect(() => {
    const safeValue = sequence.sequence || "";

    if (safeValue) {
      const validationResult = validateDnaSequence(safeValue, true);
      if (validationResult.isValid) {
        setValidation({ message: "", isAdvisory: false });
      } else {
        setValidation({
          message: validationResult.message,
          isAdvisory: validationResult.isAdvisory,
        });
      }
    } else {
      setValidation({ message: "", isAdvisory: false });
    }
  }, [sequence.sequence, index]);

  const handleSequenceChange = (e) => {
    updateSequence("sequence", e.target.value);
  };

  const handlePrimerNameChange = (e) => {
    updateSequence("primerName", e.target.value);
  };

  const handleMtkPartChange = (partType) => (e) => {
    const newValue = e.target.value;
    if (useSeparateParts) {
      updateSequence(partType, newValue);
    } else {
      updateSequence("mtkPartLeft", newValue);
      updateSequence("mtkPartRight", newValue);
    }
  };

  const handleToggleChange = (e) => {
    setUseSeparateParts(e.target.checked);
  };

  return (
    <div className="sequence-tab-content">
      <div className="form-row primer-name-row">
        <label htmlFor={`primer-name-${index}`} className="primer-name-label">
          Primer Name:
        </label>
        <input
          type="text"
          id={`primer-name-${index}`}
          value={sequence.primerName}
          onChange={handlePrimerNameChange}
          className="form-control primer-name-input"
          placeholder="Enter primer name"
        />
      </div>

      <div className="form-row sequence-header">
        <label htmlFor={`sequence-${index}`}>Sequence {index + 1}:</label>
        {charCount > 0 && (
          <div className="char-count">Length: {charCount} bp</div>
        )}
      </div>
      <div className="form-row">
        <SequenceInput
          id={`sequence-${index}`}
          value={sequence.sequence}
          onChange={handleSequenceChange}
          placeholder="Paste your DNA sequence here"
        />
        {validation.message && (
          <div
            className={`validation-message ${
              validation.isAdvisory ? "advisory" : "error"
            }`}
          >
            {validation.message}
          </div>
        )}
      </div>

      <div className="form-row mtk-parts-row">
        <div className="mtk-part-container">
          <label htmlFor={`mtk-part-left-${index}`} className="mtk-part-label">
            {useSeparateParts ? "MTK Part Number Left:" : "MTK Part Number:"}
          </label>
          <select
            id={`mtk-part-left-${index}`}
            value={sequence.mtkPartLeft || ""}
            onChange={handleMtkPartChange("mtkPartLeft")}
            className="form-control mtk-part-select"
          >
            {mtkPartOptions.map((part) => (
              <option key={part} value={part}>
                {part}
              </option>
            ))}
          </select>
        </div>
        {useSeparateParts && (
          <div className="mtk-part-container">
            <label
              htmlFor={`mtk-part-right-${index}`}
              className="mtk-part-label"
            >
              MTK Part Number Right:
            </label>
            <select
              id={`mtk-part-right-${index}`}
              value={sequence.mtkPartRight || ""}
              onChange={handleMtkPartChange("mtkPartRight")}
              className="form-control mtk-part-select"
            >
              {mtkPartOptions.map((part) => (
                <option key={part} value={part}>
                  {part}
                </option>
              ))}
            </select>
          </div>
        )}
      </div>

      <div className="form-row toggle-container">
        <label className="toggle-label">
          <input
            type="checkbox"
            checked={useSeparateParts}
            onChange={handleToggleChange}
          />
          Use separate left/right part numbers
        </label>
      </div>
    </div>
  );
}

export default SequenceTab;
