import React, { useRef, useCallback } from "react";
import "../../styles/Settings.css";

const Settings = ({
  show,
  onClose,
  formData,
  updateFormData,
  availableSpecies = [],
  settingsToggleRef,
}) => {
  const modalContentRef = useRef(null);

  const getSettingsPosition = useCallback(() => {
    if (!settingsToggleRef?.current) return {};

    const toggleRect = settingsToggleRef.current.getBoundingClientRect();
    const sidebarRect = settingsToggleRef.current
      .closest(".sidebar")
      ?.getBoundingClientRect();

    if (!sidebarRect) return {};

    return {
      position: "absolute",
      top: `${toggleRect.bottom - sidebarRect.top}px`, // Position under the button
      left: `${sidebarRect.left}px`, // Align with toggle button
      width: "90%",
      zIndex: 10,
    };
  }, [settingsToggleRef]);

  const handleOutsideClick = useCallback(
    (e) => {
      if (
        modalContentRef.current &&
        !modalContentRef.current.contains(e.target)
      ) {
        onClose();
      }
    },
    [onClose]
  );

  if (!show) return null;

  const speciesValue =
    formData.species ||
    (availableSpecies.length > 0 ? availableSpecies[0] : "");

  return (
    <div
      className="settings-modal"
      style={getSettingsPosition()}
      onClick={handleOutsideClick}
    >
      <div className="settings-content" ref={modalContentRef}>
        {/* Species Selection */}
        <div className="form-group">
          <label htmlFor="species-select">Species:</label>
          <select
            id="species-select"
            className="form-control"
            value={speciesValue}
            onChange={(e) =>
              updateFormData({ ...formData, species: e.target.value })
            }
          >
            {availableSpecies.length > 0 ? (
              availableSpecies.map((species) => (
                <option key={species} value={species}>
                  {species}
                </option>
              ))
            ) : (
              <option value="" disabled>
                No species available
              </option>
            )}
          </select>
        </div>

        {/* Kozak Selection */}
        <div className="form-group">
          <label htmlFor="kozak-select">Kozak:</label>
          <select
            id="kozak-select"
            className="form-control"
            value={formData.kozak || ""}
            onChange={(e) =>
              updateFormData({ ...formData, kozak: e.target.value })
            }
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
            value={formData.max_mut_per_site || 1}
            onChange={(e) =>
              updateFormData({
                ...formData,
                max_mut_per_site: parseInt(e.target.value, 10),
              })
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
            value={formData.results_limit || 1}
            onChange={(e) =>
              updateFormData({
                ...formData,
                results_limit:
                  e.target.value === "all"
                    ? "all"
                    : parseInt(e.target.value, 10),
              })
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
            checked={!!formData.verbose_mode}
            onChange={(e) =>
              updateFormData({ ...formData, verbose_mode: e.target.checked })
            }
          />
          <label className="form-check-label" htmlFor="verbose-mode">
            Verbose
          </label>
        </div>
      </div>
    </div>
  );
};

export default Settings;
