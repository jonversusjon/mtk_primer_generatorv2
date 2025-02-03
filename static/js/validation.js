// validation.js

document.addEventListener('DOMContentLoaded', function () {
    const form = document.getElementById('primer-form');

    form.addEventListener('input', function (event) {
        const element = event.target;

        if (element.id === 'templateSequence') {
            validateTemplateSequence(element);
        } else if (element.classList.contains('dynamic-sequence-input')) {
            validateDynamicSequence(element);
        } else if (element.id === 'numSequences') {
            // Example: Validate numSequences is not greater than 10
            if (parseInt(element.value) > 10) {
              displayValidationMessage(element, "Maximum 10 sequences allowed", false);
            } else {
              displayValidationMessage(element, "", true);
            }
        }
    });

    function validateTemplateSequence(element) {
        // Regular expression to allow for degenerate bases
        const isValid = /^[ATGCWSMKRYBDHVN]+$/i.test(element.value);

        if (isValid) {
            displayValidationMessage(element, "", true); // Clear the message if valid
        } else {
            displayValidationMessage(element, "Invalid DNA sequence. Only standard and degenerate bases (A, T, G, C, W, S, M, K, R, Y, B, D, H, V, N) are allowed.", false);
        }
    }

    function displayValidationMessage(element, message, isValid) {
        // Get the parent node of the input element
        const parentNode = element.parentNode;

        // Try to find an existing validation message
        let validationMessage = parentNode.querySelector('.validation-message');

        // If a validation message exists
        if (validationMessage) {
            // If the field is valid, remove the validation message
            if (isValid) {
                parentNode.removeChild(validationMessage);
            } else {
                // Update the existing message text
                validationMessage.textContent = message;
            }
        } else if (!isValid) {
            // Create and display a new validation message
            validationMessage = document.createElement('div');
            validationMessage.classList.add('validation-message');
            validationMessage.textContent = message;
            validationMessage.style.color = 'red'; // You can still style it as needed
            parentNode.appendChild(validationMessage);
        }
    }
    
    function validateDynamicSequence(element) {
        // Example: Check that the sequence is not empty
        if (element.value.trim() === "") {
            displayValidationMessage(element, "Sequence cannot be empty", false);
        } else {
            displayValidationMessage(element, "", true);
        }
    }

    function displayValidationMessage(element, message, isValid) {
        // Remove any existing validation message
        let validationMessage = element.parentNode.querySelector('.validation-message');
        if (validationMessage) {
          validationMessage.remove();
        }
    
        // If the field is invalid and there's a message, display it
        if (!isValid && message) {
          validationMessage = document.createElement('div');
          validationMessage.classList.add('validation-message');
          validationMessage.textContent = message;
          validationMessage.style.color = isValid ? 'green' : 'red';
          element.parentNode.appendChild(validationMessage);
        }
      }
    
    // You might also want to handle the change event for file inputs
    const fileInput = document.getElementById('fileUpload');
    if (fileInput) {
        fileInput.addEventListener('change', function (event) {
            // Add file validation logic here
            // ...
        });
    }
});