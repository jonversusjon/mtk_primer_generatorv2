import React from "react";
import "../../styles/Settings.css";

function Settings({
  show,
  onClose,
  formData,
  updateFormData,
  availableSpecies,
}) {
  if (!show) return null;

  return (
    <div className="settings-modal">
      <div className="settings-content">
        {/* Species Selection */}
        <div className="form-group">
          <label htmlFor="species-select">Species:</label>
          <select
            id="species-select"
            className="form-control"
            value={formData.species}
            onChange={(e) => updateFormData("species", e.target.value)}
          >
            <option value="">Select species...</option>
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

        <div className="settings-footer">
          <button type="button" className="btn btn-primary" onClick={onClose}>
            Close Settings
          </button>
        </div>
      </div>
    </div>
  );
}

export default Settings;
