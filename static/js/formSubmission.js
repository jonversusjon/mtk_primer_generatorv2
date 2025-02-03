// In static/js/formSubmission.js

function collectFormData() {
    const form = document.getElementById("primer-form");
    const formData = new FormData(form);
    const packagedData = {
      numSequences: parseInt(formData.get('numSequences')),
      kozak: formData.get('kozak'),
      species: formData.get('species'),
      verbose: formData.has('verbose_mode'), // Checkbox
      templateSequence: formData.get('templateSequence'),
      sequences: [],
    };
  
    for (let i = 0; i < packagedData.numSequences; i++) {
      packagedData.sequences.push({
        primerName: formData.get(`sequences[${i}][primerName]`),
        mtkPart: formData.get(`sequences[${i}][mtkPart]`),
        sequence: formData.get(`sequences[${i}][sequence]`),
      });
    }
  
    return packagedData;
  }
  
  document.addEventListener("DOMContentLoaded", function () {
    const form = document.getElementById("primer-form");
  
    form.addEventListener("submit", function (event) {
      event.preventDefault(); // Prevent default form submission
  
      const packagedData = collectFormData();
      console.log('Packaged Data:', packagedData);
  
      // Trigger HTMX post with JSON data
      htmx.trigger("#primer-form", "htmx:configRequest", {
        request: {
          body: packagedData,
          method: "POST",
        }
      });
    });
  });