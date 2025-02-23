// static/utils/utils.js

/**
 * Displays or removes a validation message next to a form element.
 * @param {HTMLElement} element - The form input element.
 * @param {string} message - The validation message to display.
 * @param {boolean} isValid - Whether the input is valid (true = remove message).
 */
function displayValidationMessage(element, message, isValid) {
  if (!element || !element.parentNode) return;

  let validationMessage = element.parentNode.querySelector(".validation-message");

  if (isValid) {
      if (validationMessage) validationMessage.remove(); // Remove existing message
  } else {
      if (!validationMessage) {
          validationMessage = document.createElement("div");
          validationMessage.classList.add("validation-message", "invalid-feedback");
          element.parentNode.appendChild(validationMessage);
      }
      validationMessage.textContent = message;
  }
}

/**
* Extracts and formats validation data for an input field.
* @param {HTMLElement} element - The form input element.
* @returns {Object} Formatted validation data.
*/
function formatValidationData(element) {
  if (!element) return {};

  return {
      field: element.getAttribute("name") || "",
      value: element.value || "",
      sequenceIndex: element.getAttribute("data-sequence-index") || null
  };
}

// Export functions for module use (if using ES6 modules)
export { displayValidationMessage, formatValidationData };
