export default class CharacterCountManager {
    constructor() {
        this.setupCharacterCounters();
    }

    /**
     * Initializes character counters for all sequence inputs.
     */
    setupCharacterCounters() {
        // Select both dynamic sequences and the template sequence
        document.querySelectorAll(".dynamic-sequence-input, .sequence-input").forEach(input => {
            input.addEventListener("input", () => this.updateCharCount(input));
            // Update character count on page load
            this.updateCharCount(input);
        });
    }

    /**
     * Updates the character count display next to the input field.
     */
    updateCharCount(input) {
        let charCountElement = null;

        // For the template sequence, which has an id "templateSequence"
        if (input.id === "templateSequence") {
            charCountElement = document.getElementById("charCount");
        } else {
            // For dynamic sequences, traverse up to find the container with class "sequence-tab-content"
            let container = input.parentElement;
            while (container && !container.classList.contains("sequence-tab-content")) {
                container = container.parentElement;
            }
            if (!container) {
                console.warn("⚠️ Could not find sequence-tab-content for:", input);
                return;
            }
            charCountElement = container.querySelector(".char-count-label");
        }

        if (charCountElement) {
            charCountElement.textContent = `Length: ${input.value.length} bp`;
            charCountElement.style.display = input.value.length ? "inline" : "none";
        } else {
            console.warn("⚠️ Character count element not found for:", input);
        }
    }
}
