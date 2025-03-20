import React, { useRef } from "react";
import "../styles/Settings.css";

const sliderValues = ["one", "a few", "many", "most", "all"];

function Settings({
  show,
  onClose,
  formData,
  updateFormData,
  availableSpecies,
}) {
  const modalContentRef = useRef(null);

  const handleOutsideClick = (e) => {
    if (
      modalContentRef.current &&
      !modalContentRef.current.contains(e.target)
    ) {
      onClose();
    }
  };

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

        {/* Results Limit Setting (Replaced Select with Slider) */}
        <div className="form-group">
          <label htmlFor="results-slider">Number of Results:</label>
          <input
            type="range"
            id="results-slider"
            className="form-control"
            min="0"
            max="4"
            step="1"
            value={sliderValues.indexOf(formData.results_limit)}
            onChange={(e) =>
              updateFormData(
                "results_limit",
                sliderValues[parseInt(e.target.value)]
              )
            }
          />
          <span style={{ marginLeft: "10px" }}>{formData.results_limit}</span>
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
      </div>
    </div>
  );
}

export default Settings;
