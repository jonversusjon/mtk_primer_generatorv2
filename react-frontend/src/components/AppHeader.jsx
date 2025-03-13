import React from "react";
import "../styles/global/theme.css";
import "../styles/AppHeader.css";

function AppHeader({ darkMode, toggleDarkMode }) {
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
        {/* Settings toggle moved to Form component */}
      </div>
    </header>
  );
}

export default AppHeader;
