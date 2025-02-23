// static/js/handlers/form-reset-handler.js

/**
 * Manages form clearing and reset functionality
 */
export default class FormResetHandler {
  constructor() {
      this.clearFormButton = document.getElementById("clearForm");

      if (this.clearFormButton) {
          this.clearFormButton.addEventListener("click", () => this.clearForm());
      }
  }

  /**
   * Clears form inputs, resets sequence handler, and clears results
   */
  clearForm() {
      // Clear file inputs safely
      const fileInputs = document.querySelectorAll("input[type='file']");
      fileInputs.forEach(input => input.value = "");

      // Reset sequenceHandlerInstance if it exists
      if (window.sequenceHandlerInstance) {
          window.sequenceHandlerInstance.currentCount = 1;
          if (typeof APP_STATE !== "undefined") {
              window.sequenceHandlerInstance.updateSequenceInputs(1, APP_STATE.testingMode);
          } else {
              console.warn("APP_STATE is not defined. Skipping updateSequenceInputs.");
          }
      }

      // Clear results section
      const resultsElement = document.getElementById("results");
      if (resultsElement) {
          resultsElement.innerHTML = "";
      }

      console.log("Form has been cleared.");
  }
}

// Initialize when DOM is loaded
document.addEventListener("DOMContentLoaded", () => {
  new FormResetHandler();
});
