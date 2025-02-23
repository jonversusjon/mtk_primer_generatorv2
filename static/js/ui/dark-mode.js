// static/js/ui/dark-mode.js

/**
 * Applies dark mode based on the user's stored preference.
 */
const applyDarkModeFromStorage = () => {
  const isDarkModeEnabled = localStorage.getItem("dark-mode") === "enabled";
  document.body.classList.toggle("dark-mode", isDarkModeEnabled);
};

/**
 * Toggles dark mode and updates the preference in local storage.
 */
const toggleDarkMode = () => {
  const isDarkModeNowEnabled = document.body.classList.toggle("dark-mode");
  localStorage.setItem("dark-mode", isDarkModeNowEnabled ? "enabled" : "disabled");
  console.log(`Dark mode ${isDarkModeNowEnabled ? "enabled" : "disabled"}.`);
};

/**
 * Initializes dark mode functionality:
 * - Applies the stored preference on page load
 * - Attaches event listener to the dark mode toggle button (if present)
 */
const initDarkMode = () => {
  applyDarkModeFromStorage();

  const darkModeButton = document.getElementById("darkModeToggle");
  if (darkModeButton) {
    darkModeButton.addEventListener("click", toggleDarkMode);
  }
};

// Initialize after DOM is fully loaded
document.addEventListener("DOMContentLoaded", initDarkMode);
