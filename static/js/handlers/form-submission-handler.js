// static/js/handlers/form-submission-handler.js

/**
 * Manages passing form data to the Flask backend using HTMX
 */
class FormSubmissionHandler {
    constructor() {
        this.form = document.getElementById("primer-form");
        this.submitButton = document.getElementById("runDesignPrimerBtn");

        if (this.form) {
            this.setupEventListeners();
        }
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
     * Updates the submit button's state based on validation messages
     */
    updateSubmitButton() {
        const errorMessages = this.form.querySelectorAll('.error-message');
        const hasErrors = Array.from(errorMessages).some(el => el.textContent.trim() !== '');
        this.submitButton.disabled = hasErrors;
    }

    /**
     * Sets up HTMX-related event listeners for form submission
     */
    setupEventListeners() {
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