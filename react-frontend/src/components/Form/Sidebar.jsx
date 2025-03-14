import React from "react";
import { Tooltip } from "react-tooltip";
import "../../styles/global/theme.css";
import "../../styles/Sidebar.css";
import Settings from "./Settings";

const SidebarItem = ({ seq, index, errors }) => (
  <div className="sidebar-item">
    {errors.length > 0 ? (
      <>
        <div
          data-tooltip-id={`error-tooltip-${index}`}
          data-tooltip-content={errors.join("\n")}
          className="error-badge-container"
        >
          <span className="error-badge">{errors.length}</span>
        </div>
        <Tooltip
          id={`error-tooltip-${index}`}
          place="right"
          style={{ whiteSpace: "pre-line" }}
        />
      </>
    ) : (
      <div className="checkmark-container">✔</div>
    )}

    <div className="sidebar-item-info">
      <div className="primer-name">{seq.primerName || "Unnamed Primer"}</div>
      <div className="sequence-length">{(seq.sequence || "").length} bp</div>
    </div>
  </div>
);

const SettingsToggle = ({
  settingsToggleRef,
  showSettings,
  setShowSettings,
}) => (
  <button
    ref={settingsToggleRef}
    type="button"
    className="settings-toggle"
    onClick={() => setShowSettings((prev) => !prev)}
    aria-label={showSettings ? "Close settings" : "Open settings"}
  >
    <span className="icon">⚙️</span>
    <span className="button-text">Settings</span>
  </button>
);

const Sidebar = ({
  sequences,
  errorsBySequence,
  settingsToggleRef,
  setShowSettings,
  showSettings,
  formData,
  setFormData,
}) => (
  <div className="sidebar">
    {/* Settings Component */}
    <Settings
      show={showSettings}
      onClose={() => setShowSettings(false)}
      formData={formData}
      updateFormData={setFormData}
      availableSpecies={formData.availableSpecies || []}
      settingsToggleRef={settingsToggleRef}
    />

    {/* Settings Toggle Button */}
    <SettingsToggle
      settingsToggleRef={settingsToggleRef}
      showSettings={showSettings}
      setShowSettings={setShowSettings}
    />

    <div className="sidebar-header">
      <h3>Sequences Overview</h3>
      <ul>
        {sequences.map((seq, index) => (
          <SidebarItem
            key={index}
            seq={seq}
            index={index}
            errors={errorsBySequence[index] || []}
          />
        ))}
      </ul>
    </div>
  </div>
);

export default Sidebar;
