// Import application state
import "./utils/app-state.js";

// Import form handlers
import "./handlers/form-reset-handler.js";
import "./handlers/form-submission-handler.js";

// Import UI components
import "./ui/dark-mode.js";
import "./ui/character-count.js";
import "./ui/sequence-tabs.js";
import "./ui/tooltip.js";

// Import utilities
import "./utils/utils.js";
import "./utils/validation.js";

import { initializeTabs } from "./ui/sequence-tabs.js";

document.addEventListener("DOMContentLoaded", () => {
    console.log("✅ Main script loaded!");

    // Initialize core components
    const sequenceHandlerInstance = new SequenceHandler();
    window.sequenceHandlerInstance = sequenceHandlerInstance;

    new FormSubmissionHandler();
    new TooltipManager();

    // Ensure character counters update after form fields are pre-filled
    setTimeout(() => {
        document.querySelectorAll(".dynamic-sequence-input").forEach(input => {
            if (input.value) {
                sequenceHandlerInstance.updateCharCount(input);
            }
        });
    }, 50);

    // Enable verbose mode if testing
    if (typeof APP_STATE !== "undefined" && APP_STATE.testingMode) {
        const verboseModeCheckbox = document.getElementById("verbose_mode");
        if (verboseModeCheckbox) {
            verboseModeCheckbox.checked = true;
        }
    }

    // Initialize tab functionality
    initializeTabs();

    // Sequence count display and controls
    const numDisplay = document.getElementById("numDisplay");
    const incrementBtn = document.getElementById("incrementBtn");
    const decrementBtn = document.getElementById("decrementBtn");

    if (incrementBtn) {
        incrementBtn.addEventListener("click", () => {
            sequenceHandlerInstance.incrementSequenceCount();
            if (numDisplay) {
                numDisplay.textContent = sequenceHandlerInstance.currentCount;
            }
        });
    }
    if (decrementBtn) {
        decrementBtn.addEventListener("click", () => {
            sequenceHandlerInstance.decrementSequenceCount();
            if (numDisplay) {
                numDisplay.textContent = sequenceHandlerInstance.currentCount;
            }
        });
    }
});
