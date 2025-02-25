// static/js/handlers/form-submission-handler.js

/**
 * Manages passing form submission and validation state
 */
export default class FormSubmissionHandler {
    constructor() {
        this.form = document.getElementById("primer-form");
        this.submitButton = document.getElementById("runDesignPrimerBtn");
        this.loadingIndicator = document.getElementById("loading-indicator");
        this.formTouched = false;
        
        if (this.form && this.submitButton) {
            this.setupFormValidation();
        }
        
        this.disableSubmitButton("Fill in required fields");
    }
    
    
    /**
     * Disables the submit button and adds a tooltip/title
     */
    disableSubmitButton(reason = "Processing...") {
        this.submitButton.disabled = true;
        this.submitButton.title = reason;
    }
    
    /**
     * Enables the submit button
     */
    enableSubmitButton() {
        this.submitButton.disabled = false;
        this.submitButton.title = "Submit form";
    }
    
    /**
     * Updates the submit button based on form validity
     */
    updateSubmitButtonState() {
        if (!this.formTouched) {
            this.disableSubmitButton("Form unchanged");
            return;
        }

        // Check if form is valid
        const formIsValid = this.checkFormValidity();

        if (formIsValid) {
            this.enableSubmitButton();
            console.log("Validation success!");

        } else {
            console.log("Validation fail: the form is not fully validated yet.");
            this.disableSubmitButton("Please correct errors");
        }
    }
    
    /**
     * Checks overall form validity by examining all validation error messages and required fields
     */
    checkFormValidity() {
        // Check for error messages
        const errorMessages = this.form.querySelectorAll('.error-message');
        const hasErrors = Array.from(errorMessages).some(el => el.textContent.trim() !== '');

        if (hasErrors) {
            return false;
        }

        // Get current number of sequences
        const numSequences = parseInt(document.getElementById('numSequencesInput').value) || 0;

        // Check required fields for each active sequence
        for (let i = 0; i < numSequences; i++) {
            const sequenceTab = this.form.querySelector(`.sequence-tab-content[data-tab-index="${i}"]`);
            if (!sequenceTab || sequenceTab.classList.contains('hidden')) continue;

            const sequence = sequenceTab.querySelector(`textarea[name="sequences[${i}][sequence]"]`);
            const primerName = sequenceTab.querySelector(`input[name="sequences[${i}][primerName]"]`);
            const mtkPart = sequenceTab.querySelector(`select[name="sequences[${i}][mtkPart]"]`);

            if (!sequence || !sequence.value.trim() ||
                !mtkPart || !mtkPart.value.trim()) {
                return false;
            }
        }

        return true;
    }
    
    /**
     * Sets up the initial form validation state
     */
    setupFormValidation() {
        // Mark form as touched when any input changes
        const inputs = this.form.querySelectorAll('input, textarea, select');
        inputs.forEach(input => {
            input.addEventListener('input', () => {
                this.formTouched = true;
                this.updateSubmitButtonState();
            });

            // Also mark as touched on focus (first interaction)
            input.addEventListener('focus', () => {
                if (!this.formTouched) {
                    this.formTouched = true;
                    this.updateSubmitButtonState();
                }
            });
        });
    }
    
    /**
     * Formats validation data for HTMX requests
     */
    formatValidationData(elt) {
        return {
            fieldName: elt.getAttribute('name'),
            fieldValue: elt.value
        };
    }
}