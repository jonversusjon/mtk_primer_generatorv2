// In a new file, e.g., static/js/utils.js
function displayValidationMessage(element, message, isValid) {
    const parentNode = element.parentNode;
    let validationMessage = parentNode.querySelector(".validation-message");
  
    if (validationMessage) {
      if (isValid) {
        parentNode.removeChild(validationMessage);
      } else {
        validationMessage.textContent = message;
        validationMessage.classList.add("invalid-feedback"); // Add a class for styling
      }
    } else if (!isValid) {
      validationMessage = document.createElement("div");
      validationMessage.classList.add("validation-message", "invalid-feedback");
      validationMessage.textContent = message;
      parentNode.appendChild(validationMessage);
    }
  }