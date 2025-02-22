// dark-mode.js

// Helper function to apply dark mode based on stored preference.
function applyDarkModeFromStorage() {
  if (localStorage.getItem("dark-mode") === "enabled") {
    document.body.classList.add("dark-mode");
  } else {
    document.body.classList.remove("dark-mode");
  }
}

//  function to toggle dark mode and update local storage.
function toggleDarkMode() {
  document.body.classList.toggle("dark-mode");
  const isDarkMode = document.body.classList.contains("dark-mode");
  localStorage.setItem("dark-mode", isDarkMode ? "enabled" : "disabled");
}

// Initialize dark mode functionality:
// - Apply stored preference
// - Attach event listener to the dark mode toggle button (if it exists)
function initDarkMode() {
  applyDarkModeFromStorage();

  const darkModeButton = document.getElementById("darkModeToggle");
  if (darkModeButton) {
    darkModeButton.addEventListener("click", toggleDarkMode);
  }
}

// Run the initialization once the DOM content has loaded.
document.addEventListener("DOMContentLoaded", initDarkMode);