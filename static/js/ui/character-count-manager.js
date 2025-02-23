export default class CharacterCountManager {
    constructor() {
        this.setupCharacterCounters();
    }

    /**
     * Initializes character counters for all sequence inputs.
     */
    setupCharacterCounters() {
        document.querySelectorAll(".dynamic-sequence-input").forEach(input => {
            input.addEventListener("input", () => this.updateCharCount(input));
            // Ensure character count is updated on page load
            this.updateCharCount(input);
        });
    }

    /**
     * Updates the character count display next to the input field.
     */
    updateCharCount(input) {
        // Traverse up the DOM to find the sequence-tab-content wrapper
        let container = input.parentElement;
        while (container && !container.classList.contains("sequence-tab-content")) {
            container = container.parentElement;
        }
    
        if (!container) {
            console.warn("⚠️ Could not find sequence-tab-content for:", input);
            return;
        }
    
        // Now, look for the character count element inside this container
        const charCountElement = container.querySelector(".char-count-label");
    
        if (charCountElement) {
            charCountElement.textContent = `Length: ${input.value.length} bp`;
            charCountElement.style.display = input.value.length ? "inline" : "none";
        } else {
            console.warn("⚠️ Character count element not found for:", input);
        }
    }
    
}