import React from "react";

function AppHeader({ darkMode, toggleDarkMode }) {
  return (
    <div className="app-header">
      <h1>MTK Primer Domesticator</h1>

      <div className="header-controls">
        <button
          className="dark-mode-btn"
          onClick={toggleDarkMode}
          aria-label={darkMode ? "Switch to light mode" : "Switch to dark mode"}
        >
          <span className="moon-icon">{darkMode ? "☀️" : "🌙"}</span>
          <span className="button-text">
            {darkMode ? "Light Mode" : "Dark Mode"}
          </span>
        </button>
      </div>
    </div>
  );
}

export default AppHeader;
