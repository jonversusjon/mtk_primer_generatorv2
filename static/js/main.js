// static/js/main.js

// Import application state
import "./utils/app-state.js";

// Import form handlers
import FormResetHandler from "./handlers/form-reset-handler.js";
import FormSubmissionHandler from "./handlers/form-submission-handler.js";

// Import UI components
import DarkModeManager from "./ui/dark-mode-manager.js";
import CharacterCountManager from "./ui/character-count-manager.js";
import TabScroller from "./ui/sequence-tabs/sequence-tabs-scroller.js";
import SequenceTabsManager from "./ui/sequence-tabs/sequence-tabs-manager.js";
import TooltipManager from "./ui/tooltip-manager.js";

// Import utilities
import "./utils/utils.js";
import "./utils/validation.js";

document.addEventListener("DOMContentLoaded", () => {
    console.log("Frontend reporting for business!");

    // Initialize core components
    const sequenceTabsManagerInstance = new SequenceTabsManager();
    window.sequenceTabsManagerInstance = sequenceTabsManagerInstance;

    new FormSubmissionHandler();
    new FormResetHandler();
    new TabScroller("sequence-tabs-container");
    new CharacterCountManager();
    new TooltipManager();
    new DarkModeManager();

    // Ensure character counters update after form fields are pre-filled
    setTimeout(() => {
        document.querySelectorAll(".dynamic-sequence-input").forEach(input => {
            if (input.value && sequenceTabsManagerInstance.updateCharCount) {
                sequenceTabsManagerInstance.updateCharCount(input);
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

    // Sequence count display and controls
    const numDisplay = document.getElementById("numDisplay");
    const incrementBtn = document.getElementById("incrementBtn");
    const decrementBtn = document.getElementById("decrementBtn");

    if (incrementBtn) {
        incrementBtn.addEventListener("click", () => {
            sequenceTabsManagerInstance.incrementSequenceCount();
            if (numDisplay) {
                numDisplay.textContent = sequenceTabsManagerInstance.currentCount;
            }
        });
    }
    if (decrementBtn) {
        decrementBtn.addEventListener("click", () => {
            sequenceTabsManagerInstance.decrementSequenceCount();
            if (numDisplay) {
                numDisplay.textContent = sequenceTabsManagerInstance.currentCount;
            }
        });
    }
});
