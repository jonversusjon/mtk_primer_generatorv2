// dark-mode.js
export function toggleDarkMode() {
    document.body.classList.toggle("dark-mode");
    let isDarkMode = document.body.classList.contains("dark-mode");
    localStorage.setItem("dark-mode", isDarkMode ? "enabled" : "disabled");
}

document.addEventListener("DOMContentLoaded", function () {
    const darkModeButton = document.getElementById("darkModeToggle");
    if (darkModeButton) {
        darkModeButton.addEventListener("click", toggleDarkMode);
    }

    if (localStorage.getItem("dark-mode") === "enabled") {
        document.body.classList.add("dark-mode");
    }
});

