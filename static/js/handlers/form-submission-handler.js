// static/js/handlers/form-submission-handler.js

/**
 * Manages passing form submission and and validation state
 */
export default class FormSubmissionHandler {
    constructor() {
        this.form = document.getElementById("primer-form");
        this.submitButton = document.getElementById("runDesignPrimerBtn");
        this.loadingIndicator = document.getElementById("loading-indicator");
        this.formTouched = false;


        if (this.form && this.submitButton) {
            this.setupEventListeners();
            this.setupFormValidation();
        }

        this.disableSubmitButton("Fill in required fields");
    }

    /**
     * Collects and structures form data for submission
     */
    collectFormData() {
        const formData = new FormData(this.form);
        const verboseModeCheckbox = this.form.querySelector('input[name="verbose_mode"]');

        const data = {
            numSequences: parseInt(formData.get('numSequences')) || 0,
            kozak: formData.get('kozak') || '',
            species: formData.get('species') || '',
            'max-mut-per-site': formData.get('max-mut-per-site') || '3',
            verboseMode: verboseModeCheckbox ? verboseModeCheckbox.checked : false,
            templateSequence: formData.get('templateSequence') || '',
            sequences: []
        };

        // Build sequences array
        for (let i = 0; i < data.numSequences; i++) {
            data.sequences.push({
                primerName: formData.get(`sequences[${i}][primerName]`) || '',
                mtkPart: formData.get(`sequences[${i}][mtkPart]`) || '',
                sequence: formData.get(`sequences[${i}][sequence]`) || ''
            });
        }

        console.log('Collected form data:', data);
        return data;
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
        } else {
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
     * Sets up HTMX-related event listeners for form submission
     */
    setupEventListeners() {
        // Handle form submission
        this.form.addEventListener('submit', (event) => {
            // Prevent default form submission if already processing
            if (this.submitButton.disabled && this.submitButton.classList.contains('processing')) {
                event.preventDefault();
                event.stopPropagation();
                return false;
            }
            
            // Disable button and show loading state
            this.disableSubmitButton();
            this.submitButton.classList.add('processing');
            this.submitButton.innerHTML = '<span class="spinner"></span> Processing...';
        });

        // Handle HTMX request configuration
        document.body.addEventListener('htmx:configRequest', (evt) => {
            const elt = evt.detail.elt;
            const isValidationRequest = elt.getAttribute('hx-post')?.includes('validate_field');

            console.log('HTMX Request Config:', {
                element: elt,
                isValidation: isValidationRequest,
                currentHeaders: evt.detail.headers,
                currentParams: evt.detail.parameters
            });

            if (isValidationRequest) {
                evt.detail.parameters = this.formatValidationData(elt);
            } else if (elt.id === 'primer-form') {
                evt.detail.parameters = this.collectFormData();
            }
        });

        // Listen for validation responses and update the submit button
        document.body.addEventListener('htmx:afterOnLoad', (evt) => {
            if (evt.detail.elt.getAttribute('hx-post')?.includes('validate_field')) {
                this.updateSubmitButton();
            }
        });

        // Logging before request
        document.body.addEventListener('htmx:beforeRequest', (evt) => {
            console.log('HTMX Request:', evt.detail.parameters);
        });

        // Logging after request
        document.body.addEventListener('htmx:afterRequest', (evt) => {
            console.log('HTMX Response:', {
                status: evt.detail.xhr.status,
                response: evt.detail.xhr.response,
                headers: evt.detail.xhr.getAllResponseHeaders()
            });
        });

        // Error handling for HTMX requests
        document.body.addEventListener('htmx:responseError', (evt) => {
            console.error('HTMX Response Error:', {
                error: evt.detail.error,
                xhr: evt.detail.xhr,
                status: evt.detail.xhr.status,
                response: evt.detail.xhr.response
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

// Initialize the handler when the DOM is loaded
document.addEventListener("DOMContentLoaded", () => {
    new FormSubmissionHandler();
});