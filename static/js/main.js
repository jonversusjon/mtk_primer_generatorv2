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
    const sequenceTabsManager = new SequenceTabsManager();
    window.sequenceTabsManagerInstance = sequenceTabsManager;

    new FormSubmissionHandler();
    new FormResetHandler();
    new TabScroller("sequence-tabs-container");
    new CharacterCountManager();
    new TooltipManager();
    new DarkModeManager();

    // Update character counters for pre-filled fields
    setTimeout(() => {
        document.querySelectorAll(".dynamic-sequence-input").forEach(input => {
            if (input.value) {
                sequenceTabsManager.updateCharCount?.(input);
            }
        });
    }, 50);

    const verboseModeCheckbox = document.getElementById("verbose_mode");
    if (verboseModeCheckbox) {
        verboseModeCheckbox.checked = true;
    }


    // Sequence count display and controls
    const numDisplay = document.getElementById("numDisplay");
    const updateNumDisplay = () => {
        if (numDisplay) {
            numDisplay.textContent = sequenceTabsManager.currentCount;
        }
    };

    document.getElementById("incrementBtn")?.addEventListener("click", () => {
        sequenceTabsManager.incrementSequenceCount();
        updateNumDisplay();
    });

    document.getElementById("decrementBtn")?.addEventListener("click", () => {
        sequenceTabsManager.decrementSequenceCount();
        updateNumDisplay();
    });
});
