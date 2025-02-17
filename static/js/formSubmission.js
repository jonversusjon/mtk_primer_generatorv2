// Static/js/formSubmission.js

// Helper function to format field validation data
function formatValidationData(element) {
    const field = element.getAttribute('name');
    const sequenceIndex = element.getAttribute('data-sequence-index');

    return {
        field: field,
        value: element.value,
        sequenceIndex: sequenceIndex
    };
}

// Collect all form data for submission
function collectFormData() {
    const form = document.getElementById("primer-form");
    const formData = new FormData(form);
    
    // Get checkbox value directly from the element
    const verboseModeCheckbox = form.querySelector('input[name="verbose_mode"]');
    
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
    const numSeq = data.numSequences;
    for (let i = 0; i < numSeq; i++) {
        data.sequences.push({
            primerName: formData.get(`sequences[${i}][primerName]`) || '',
            mtkPart: formData.get(`sequences[${i}][mtkPart]`) || '',
            sequence: formData.get(`sequences[${i}][sequence]`) || ''
        });
    }
    
    console.log('Collected form data:', data);
    return data;
}

// Handle HTMX requests configuration
document.body.addEventListener('htmx:configRequest', function(evt) {
    const elt = evt.detail.elt;
    const isValidationRequest = elt.getAttribute('hx-post')?.includes('validate_field');
    
    // Log the event details
    console.log('HTMX Request Config:', {
        element: elt,
        isValidation: isValidationRequest,
        currentHeaders: evt.detail.headers,
        currentParams: evt.detail.parameters
    });

    if (isValidationRequest) {
        const validationData = formatValidationData(elt);
        console.log('Validation Data:', validationData);
        evt.detail.parameters = validationData;
    } 
    else if (elt.id === 'primer-form') {
        const formData = collectFormData();
        console.log('Form Submission Data:', formData);
        evt.detail.parameters = formData;
    }
});

// Handle form submission
document.addEventListener("DOMContentLoaded", function () {
    const form = document.getElementById("primer-form");
    const submitButton = document.getElementById("runDesignPrimerBtn");

    // Enable/disable submit button based on form validity
    function updateSubmitButton() {
        const errorMessages = form.querySelectorAll('.error-message');
        const hasErrors = Array.from(errorMessages).some(el => el.textContent.trim() !== '');
        submitButton.disabled = hasErrors;
    }

    // Listen for validation responses
    document.body.addEventListener('htmx:afterOnLoad', function (evt) {
        if (evt.detail.elt.getAttribute('hx-post')?.includes('validate_field')) {
            updateSubmitButton();
        }
    });

    // Log responses for debugging
    document.body.addEventListener('htmx:beforeRequest', function (evt) {
        console.log('Sending request:', evt.detail.parameters);
    });

    document.body.addEventListener('htmx:afterRequest', function (evt) {
        console.log('Received response:', evt.detail.xhr.response);
    });
});

// Add more detailed error logging
document.body.addEventListener('htmx:beforeRequest', function(evt) {
    console.log('Before Request:', {
        url: evt.detail.path,
        method: evt.detail.method,
        headers: evt.detail.headers,
        parameters: evt.detail.parameters
    });
});

document.body.addEventListener('htmx:beforeSend', function(evt) {
    console.log('Request details before send:', {
        path: evt.detail.path,
        method: evt.detail.method,
        headers: evt.detail.headers,
        parameters: evt.detail.parameters
    });
});

document.body.addEventListener('htmx:afterRequest', function(evt) {
    console.log('After Request:', {
        status: evt.detail.xhr.status,
        response: evt.detail.xhr.response,
        headers: evt.detail.xhr.getAllResponseHeaders()
    });
});

// Add error handling
document.body.addEventListener('htmx:responseError', function(evt) {
    console.error('HTMX Response Error:', {
        error: evt.detail.error,
        xhr: evt.detail.xhr,
        status: evt.detail.xhr.status,
        response: evt.detail.xhr.response
    });
});

// Handle form submission
document.addEventListener("DOMContentLoaded", function() {
    const form = document.getElementById("primer-form");
    const submitButton = document.getElementById("runDesignPrimerBtn");
    
    // Enable/disable submit button based on form validity
    function updateSubmitButton() {
        const errorMessages = form.querySelectorAll('.error-message');
        const hasErrors = Array.from(errorMessages).some(el => el.textContent.trim() !== '');
        submitButton.disabled = hasErrors;
    }
    
    // Listen for validation responses
    document.body.addEventListener('htmx:afterOnLoad', function(evt) {
        if (evt.detail.elt.getAttribute('hx-post')?.includes('validate_field')) {
            updateSubmitButton();
        }
    });
});