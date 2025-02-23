// static/js/ui/DarkMode.js

export default class DarkModeManager {
  constructor(toggleButtonId = "darkModeToggle") {
    this.toggleButton = document.getElementById(toggleButtonId);
    this.init();
  }

  /**
   * Applies dark mode based on the user's stored preference.
   */
  applyDarkModeFromStorage() {
    const isDarkModeEnabled = localStorage.getItem("dark-mode") === "enabled";
    document.body.classList.toggle("dark-mode", isDarkModeEnabled);
  }

  /**
   * Toggles dark mode and updates the preference in local storage.
   */
  toggleDarkMode() {
    const isDarkModeNowEnabled = document.body.classList.toggle("dark-mode");
    localStorage.setItem("dark-mode", isDarkModeNowEnabled ? "enabled" : "disabled");
    console.log(`Dark mode ${isDarkModeNowEnabled ? "enabled" : "disabled"}.`);
  }

  /**
   * Initializes dark mode functionality:
   * - Applies the stored preference on page load
   * - Attaches event listener to the dark mode toggle button (if present)
   */
  init() {
    this.applyDarkModeFromStorage();

    if (this.toggleButton) {
      this.toggleButton.addEventListener("click", () => this.toggleDarkMode());
    }
  }
}
