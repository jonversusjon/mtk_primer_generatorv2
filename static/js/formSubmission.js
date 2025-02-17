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
  
  return {
      numSequences: parseInt(formData.get('numSequences')),
      kozak: formData.get('kozak'),
      species: formData.get('species'),
      'max-mut-per-site': formData.get('max-mut-per-site'),
      verboseMode: formData.has('verbose_mode'),
      templateSequence: formData.get('templateSequence'),
      sequences: Array.from({ length: parseInt(formData.get('numSequences')) }, (_, i) => ({
          primerName: formData.get(`sequences[${i}][primerName]`),
          mtkPart: formData.get(`sequences[${i}][mtkPart]`),
          sequence: formData.get(`sequences[${i}][sequence]`)
      }))
  };
}

// Handle HTMX requests configuration
document.body.addEventListener('htmx:configRequest', function(evt) {
  const elt = evt.detail.elt;
  const isValidationRequest = elt.getAttribute('hx-post')?.includes('validate_field');
  
  // Set content type header
  evt.detail.headers = {
      ...evt.detail.headers,
      'Content-Type': 'application/json'
  };
  
  // Handle field validation requests
  if (isValidationRequest) {
      evt.detail.parameters = formatValidationData(elt);
      evt.detail.unfilteredParameters = formatValidationData(elt);
  } 
  // Handle form submission
  else if (elt.id === 'primer-form') {
      const formData = collectFormData();
      evt.detail.parameters = formData;
      evt.detail.unfilteredParameters = formData;
  }
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
  
  // Log responses for debugging
  document.body.addEventListener('htmx:beforeRequest', function(evt) {
      console.log('Sending request:', evt.detail.parameters);
  });
  
  document.body.addEventListener('htmx:afterRequest', function(evt) {
      console.log('Received response:', evt.detail.xhr.response);
  });
});