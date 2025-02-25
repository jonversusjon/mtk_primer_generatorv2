// static/js/handlers/form-reset-handler.js

/**
 * Manages form clearing and reset functionality.
 */
export default class FormResetHandler {
    constructor() {
        // Reference the clear button and the form.
        this.clearFormButton = document.getElementById("clearForm");
        this.form = document.getElementById("primer-form");

        // Attach the clear functionality if the clear button exists.
        if (this.clearFormButton) {
            this.clearFormButton.addEventListener("click", () => this.clearForm());
        }
    }

    /**
     * Clears all form inputs, resets the sequence handler, and clears the results section.
     */
    clearForm() {
        // Reset the entire form.
        if (this.form) {
            this.form.reset();
        }

        // Clear file inputs safely.
        const fileInputs = document.querySelectorAll("input[type='file']");
        fileInputs.forEach(input => input.value = "");

        // Reset the global sequence handler if it exists.
        if (window.sequenceHandlerInstance) {
            window.sequenceHandlerInstance.currentCount = 1;
            if (typeof APP_STATE !== "undefined") {
                window.sequenceHandlerInstance.updateSequenceInputs(1, APP_STATE.testingMode);
            } else {
                console.warn("APP_STATE is not defined. Skipping updateSequenceInputs.");
            }
        }

        // Clear the results section if present.
        const resultsElement = document.getElementById("results");
        if (resultsElement) {
            resultsElement.innerHTML = "";
        }

        console.log("Form has been cleared.");
    }
}
