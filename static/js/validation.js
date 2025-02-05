document.addEventListener("DOMContentLoaded", function () {
    const form = document.getElementById("primer-form");
    const runDesignPrimerBtn = document.getElementById("runDesignPrimerBtn");
    const numSequencesInput = document.getElementById("numSequences");
    const templateSeqEl = document.getElementById("templateSequence");

    // Store field validity state.
    // We’ll use each field’s id (or name if no id) as the key.
    let fieldValidity = {};

    // Helper: Get a unique key for an element (id if available, otherwise name)
    function getFieldKey(element) {
        return element.id || element.name;
    }

    // Update the submit button so that it is disabled if any tracked field is invalid.
    function updateSubmitButtonState() {
        let allValid = true;
        for (const key in fieldValidity) {
            if (!fieldValidity[key]) {
                allValid = false;
                break;
            }
        }
        runDesignPrimerBtn.disabled = !allValid;
    }

    // Show (or hide) the validation message for an element and update its validity.
    function displayValidationMessage(element, message) {
        let validationMessage = element.nextElementSibling;
        if (!validationMessage || !validationMessage.classList.contains("validation-message")) {
            validationMessage = document.createElement("div");
            validationMessage.classList.add("validation-message");
            validationMessage.style.color = "red";
            element.parentNode.insertBefore(validationMessage, element.nextSibling);
        }
        validationMessage.textContent = message;
        validationMessage.style.display = message ? "block" : "none";

        // Use a stable key to record this field’s validity.
        const key = getFieldKey(element);
        if (key) {
            if (message) {
                element.classList.add("validation-fail");
                fieldValidity[key] = false;
            } else {
                element.classList.remove("validation-fail");
                fieldValidity[key] = true;
            }
        }
        updateSubmitButtonState();
    }

    // Validate the template sequence (an optional field).
    function validateTemplateSequence() {
        const value = templateSeqEl.value.trim();
        const allowedRegex = /^[ATGCWSMKRYBDHVN]*$/i;
        if (value !== "" && !allowedRegex.test(value)) {
            displayValidationMessage(
                templateSeqEl,
                "Only valid DNA bases (A, T, G, C, W, S, M, K, R, Y, B, D, H, V, N) allowed."
            );
            return false;
        }
        displayValidationMessage(templateSeqEl, "");
        return true;
    }

    // Validate a dynamic (amplicon) sequence.
    function validateDynamicSequence(element) {
        const value = element.value.trim();
        console.log("Validating dynamic sequence:", value);
        const allowedRegex = /^[ATGCWSMKRYBDHVN]+$/i;
      
        if (value === "") {
          displayValidationMessage(element, "Sequence cannot be empty.");
          return false;
        }
        if (!allowedRegex.test(value)) {
          console.log("Invalid sequence detected:", value);
          displayValidationMessage(element, "Only valid DNA bases allowed (A, T, G, C, W, S, M, K, R, Y, B, D, H, V, N).");
          return false;
        }
        if (value.length < 80) {
          displayValidationMessage(element, "Amplicon sequence must be at least 80 bp.");
          return false;
        }
      
        // Check against template sequence
        const templateValue = templateSeqEl.value.trim();
        if (templateValue !== "" && (value.length > templateValue.length || !templateValue.includes(value))) {
          displayValidationMessage(element, "Amplicon sequence must be within the template.");
          return false;
        }
      
        displayValidationMessage(element, "");
        return true;
      }
      

    // Validate the MTK Part selection.
    function validateMtkPart(element) {
        if (element.value.trim() === "") {
            displayValidationMessage(element, "MTK Part number is required.");
            return false;
        }
        displayValidationMessage(element, "");
        return true;
    }

    // Validate the number of sequences input.
    function validateNumSequences() {
        const numSequences = parseInt(numSequencesInput.value) || 0;
        if (numSequences < 1 || numSequences > 10) {
            displayValidationMessage(numSequencesInput, "Number of sequences must be between 1 and 10.");
            return false;
        }
        displayValidationMessage(numSequencesInput, "");
        return true;
    }

    // Listen for input events and validate only the changed field (and any dependents).
    form.addEventListener("input", function (event) {
        console.log("Input event on:", event.target.id, "value:", event.target.value);
        const target = event.target;
        if (target.id === "templateSequence") {
            validateTemplateSequence();
            // Revalidate all dynamic sequence inputs as the template affects them.
            form.querySelectorAll(".dynamic-sequence-input").forEach(validateDynamicSequence);
        } else if (target.classList.contains("dynamic-sequence-input")) {
            validateDynamicSequence(target);
        } else if (target.id.startsWith("mtkPart")) {
            validateMtkPart(target);
        } else if (target.id === "numSequences") {
            validateNumSequences();
        }
    });

    // Run initial validations (and initialize fieldValidity) on page load.
    validateNumSequences();
    validateTemplateSequence();
    form.querySelectorAll(".dynamic-sequence-input").forEach(validateDynamicSequence);
    form.querySelectorAll("[id^='mtkPart']").forEach(validateMtkPart);
});
