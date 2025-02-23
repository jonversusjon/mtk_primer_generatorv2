class CharacterCountManager {
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
        // Find corresponding counter element
        const charCountElement = input.closest(".character-counter-container")?.querySelector(".char-count");

        if (charCountElement) {
            charCountElement.textContent = `${input.value.length} / ${input.maxLength}`;
            charCountElement.style.display = "inline"; // Ensure it's visible
        } else {
            console.warn("⚠️ Character count element not found for:", input);
        }
    }
}

// Initialize after DOM is loaded
document.addEventListener("DOMContentLoaded", () => {
    new CharacterCountManager();
});
