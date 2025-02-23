/**
 * Manages form clearing and reset functionality
 */
class FormHandler {
    constructor() {
      this.clearFormButton = document.getElementById("clearForm");
      if (this.clearFormButton) {
        this.clearFormButton.addEventListener("click", () => this.clearForm());
      }
    }
  
    clearForm() {
      document.querySelectorAll("#fileUpload").forEach(input => input.value = "");
      if (window.sequenceHandlerInstance) {
        window.sequenceHandlerInstance.currentCount = 1;
        window.sequenceHandlerInstance.updateSequenceInputs(1, APP_STATE.testingMode);
      }
      const resultsElement = document.getElementById("results");
      if (resultsElement) {
        resultsElement.innerHTML = "";
      }
    }
  }