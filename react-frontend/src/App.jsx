// src/App.jsx
import React, { useState } from "react";
import { BrowserRouter as Router, Routes, Route } from "react-router-dom";
import AppHeader from "./components/AppHeader";
import Banner from "./components/Banner";
import FormPage from "./pages/FormPage";
import ResultsPage from "./pages/ResultsPage";
import { useDarkMode } from "./hooks/useDarkMode";

// Styles
import "./styles/app.css";
import "./styles/base.css";
import "./styles/dark-mode.css";
import "./styles/scrollbar.css";
import "./styles/theme.css";

function App() {
  const [darkMode, toggleDarkMode] = useDarkMode();
  const [showSettings, setShowSettings] = useState(false);

  // We'll store 'results' at the top-level so both pages can share it
  const [results, setResults] = useState(null);

  // Toggle advanced settings in the Form
  const toggleSettings = () => {
    setShowSettings(!showSettings);
  };

  return (
    <Router>
      <div className={`app ${darkMode ? "dark-mode" : ""}`}>
        <AppHeader
          darkMode={darkMode}
          toggleDarkMode={toggleDarkMode}
          toggleSettings={toggleSettings}
        />
        <Banner />
        <div className="app-container">
          <Routes>
            {/* 
              1) The Form Page ("/"): handles the form, calls API, 
                 and navigates to /results on success.
            */}
            <Route
              path="/"
              element={
                <FormPage
                  showSettings={showSettings}
                  setShowSettings={setShowSettings}
                  setResults={setResults}
                />
              }
            />

            {/* 
              2) The Results Page ("/results"): displays the results,
                 or redirects back to "/" if no results found. 
            */}
            <Route
              path="/results"
              element={<ResultsPage results={results} />}
            />
          </Routes>
        </div>
      </div>
    </Router>
  );
}

export default App;
