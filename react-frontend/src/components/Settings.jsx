import React, { useRef } from "react";
import "../styles/Settings.css";

function Settings({
  show,
  onClose,
  formData,
  updateFormData,
  availableSpecies,
}) {
  // Create a ref for the modal content - hooks must be called unconditionally
  const modalContentRef = useRef(null);

  // Handle clicks outside the modal content
  const handleOutsideClick = (e) => {
    // If we click outside the modal content, close the modal
    if (
      modalContentRef.current &&
      !modalContentRef.current.contains(e.target)
    ) {
      onClose();
    }
  };

  // If not shown, return null - but AFTER calling hooks
  if (!show) return null;

  return (
    <div className="settings-modal" onClick={handleOutsideClick}>
      <div className="settings-content" ref={modalContentRef}>
        {/* Species Selection */}
        <div className="form-group">
          <label htmlFor="species-select">Species:</label>
          <select
            id="species-select"
            className="form-control"
            value={formData.species}
            onChange={(e) => updateFormData("species", e.target.value)}
          >
            {availableSpecies.map((species) => (
              <option key={species} value={species}>
                {species}
              </option>
            ))}
          </select>
        </div>

        {/* Kozak Selection */}
        <div className="form-group">
          <label htmlFor="kozak-select">Kozak:</label>
          <select
            id="kozak-select"
            className="form-control"
            value={formData.kozak}
            onChange={(e) => updateFormData("kozak", e.target.value)}
          >
            <option value="MTK">MTK</option>
            <option value="Canonical">Canonical</option>
          </select>
        </div>

        {/* Mutations Setting */}
        <div className="form-group">
          <label htmlFor="mutations-select">Max mutations per site:</label>
          <select
            id="mutations-select"
            className="form-control"
            value={formData.max_mut_per_site}
            onChange={(e) =>
              updateFormData("max_mut_per_site", parseInt(e.target.value))
            }
          >
            <option value="1">1</option>
            <option value="2">2</option>
            <option value="3">3</option>
          </select>
        </div>

        {/* Results Limit Setting */}
        <div className="form-group">
          <label htmlFor="results-select">Number of Results:</label>
          <select
            id="results-select"
            className="form-control"
            value={formData.results_limit}
            onChange={(e) =>
              updateFormData(
                "results_limit",
                e.target.value === "all" ? "all" : parseInt(e.target.value)
              )
            }
          >
            {[...Array(6).keys()].map((i) => (
              <option key={i} value={2 ** i}>
                {2 ** i}
              </option>
            ))}
            <option value="all">All</option>
          </select>
        </div>

        {/* Verbose Mode Toggle */}
        <div className="form-check">
          <input
            type="checkbox"
            className="form-check-input"
            id="verbose-mode"
            checked={formData.verbose_mode}
            onChange={(e) => updateFormData("verbose_mode", e.target.checked)}
          />
          <label className="form-check-label" htmlFor="verbose-mode">
            Verbose
          </label>
        </div>

        {/* No footer with close button - will close by clicking outside or toggling settings button */}
      </div>
    </div>
  );
}

export default Settings;