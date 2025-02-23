// static/js/utils/validation.js

document.addEventListener("DOMContentLoaded", function () {
    const form = document.getElementById("primer-form");
    if (!form) return; // Prevent execution if form is missing

    const runDesignPrimerBtn = document.getElementById("runDesignPrimerBtn");
    const templateSeqEl = document.getElementById("templateSequence");

    // Store field validity state
    let fieldValidity = {};

    /**
     * Gets a unique key for a form element (uses ID if available, otherwise name)
     */
    function getFieldKey(element) {
        return element.id || element.name;
    }

    /**
     * Updates the submit button state based on field validity.
     */
    function updateSubmitButtonState() {
        runDesignPrimerBtn.disabled = Object.values(fieldValidity).some(isValid => !isValid);
    }

    /**
     * Displays or removes a validation message for an element.
     * @param {HTMLElement} element - The form field.
     * @param {string} message - The validation message (empty string removes message).
     */
    function displayValidationMessage(element, message) {
        if (!element || !element.parentNode) return;

        let validationMessage = element.parentNode.querySelector(".validation-message");

        if (message) {
            if (!validationMessage) {
                validationMessage = document.createElement("div");
                validationMessage.classList.add("validation-message", "invalid-feedback");
                element.parentNode.insertBefore(validationMessage, element.nextSibling);
            }
            validationMessage.textContent = message;
            element.classList.add("validation-fail");
            fieldValidity[getFieldKey(element)] = false;
        } else {
            if (validationMessage) validationMessage.remove();
            element.classList.remove("validation-fail");
            fieldValidity[getFieldKey(element)] = true;
        }

        updateSubmitButtonState();
    }

    /**
     * Validates the template sequence (optional field).
     */
    function validateTemplateSequence() {
        if (!templateSeqEl) return true;

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

    /**
     * Validates a dynamic sequence input field.
     * @param {HTMLElement} element - The input field.
     */
    function validateDynamicSequence(element) {
        const value = element.value.trim();
        const allowedRegex = /^[ATGCWSMKRYBDHVN]+$/i;

        if (value === "") {
            displayValidationMessage(element, "Sequence cannot be empty.");
            return false;
        }
        if (!allowedRegex.test(value)) {
            console.warn("Invalid sequence detected:", value);
            displayValidationMessage(element, "Only valid DNA bases allowed (A, T, G, C, W, S, M, K, R, Y, B, D, H, V, N).");
            return false;
        }
        if (value.length < 80) {
            displayValidationMessage(element, "Amplicon sequence must be at least 80 bp.");
            return false;
        }

        // Ensure amplicon sequence fits within template sequence
        const templateValue = templateSeqEl?.value.trim();
        if (templateValue && (value.length > templateValue.length || !templateValue.includes(value))) {
            displayValidationMessage(element, "Amplicon sequence must be within the template.");
            return false;
        }

        displayValidationMessage(element, "");
        return true;
    }

    /**
     * Validates the MTK Part selection field.
     * @param {HTMLElement} element - The select dropdown.
     */
    function validateMtkPart(element) {
        if (element.value.trim() === "") {
            displayValidationMessage(element, "MTK Part number is required.");
            return false;
        }
        displayValidationMessage(element, "");
        return true;
    }

    /**
     * Handles real-time validation on user input.
     */
    form.addEventListener("input", function (event) {
        const target = event.target;
        console.debug("Validating input:", target.id, "value:", target.value);

        if (target.id === "templateSequence") {
            validateTemplateSequence();
            form.querySelectorAll(".dynamic-sequence-input").forEach(validateDynamicSequence);
        } else if (target.classList.contains("dynamic-sequence-input")) {
            validateDynamicSequence(target);
        } else if (target.id.startsWith("mtkPart")) {
            validateMtkPart(target);
        }
    });

    /**
     * Runs initial validations on page load.
     */
    validateTemplateSequence();
    form.querySelectorAll(".dynamic-sequence-input").forEach(validateDynamicSequence);
    form.querySelectorAll("[id^='mtkPart']").forEach(validateMtkPart);
});
