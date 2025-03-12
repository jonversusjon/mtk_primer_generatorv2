import React from "react";
import "../styles/AppHeader.css";

function AppHeader({ darkMode, toggleDarkMode, toggleSettings }) {
  return (
    <header className="app-header">
      <div className="header-left">
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
      <div className="header-center">
        <h1>MTK Advanced Primer Designer</h1>
      </div>
      <div className="header-right">
        <button
          className="settings-toggle"
          onClick={toggleSettings}
          aria-label="Open settings"
        >
          <span className="icon">⚙️</span>
          <span className="button-text">Settings</span>
        </button>
      </div>
    </header>
  );
}

export default AppHeader;
