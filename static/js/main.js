// static/js/main.js

// Initialize everything when DOM is ready
document.addEventListener("DOMContentLoaded", () => {
  const sequenceHandlerInstance = new SequenceHandler();
  window.sequenceHandlerInstance = sequenceHandlerInstance;

  new FormHandler();
  new TooltipManager();

  // Ensure character counters update after form fields are pre-filled
  setTimeout(() => {
    document.querySelectorAll(".dynamic-sequence-input").forEach(input => {
      if (input.value) {
        window.sequenceHandlerInstance.updateCharCount(input);
      }
    });
  }, 50);

  if (APP_STATE.testingMode) {
    document.getElementById("verbose_mode").checked = true;
  }

  initializeTabs();

  const numDisplay = document.getElementById("numDisplay");
  const incrementBtn = document.getElementById("incrementBtn");
  const decrementBtn = document.getElementById("decrementBtn");

  if (incrementBtn) {
    incrementBtn.addEventListener("click", () => {
      sequenceHandlerInstance.incrementSequenceCount();
      if (numDisplay) {
        numDisplay.textContent = sequenceHandlerInstance.currentCount;
      }
    });
  }
  if (decrementBtn) {
    decrementBtn.addEventListener("click", () => {
      sequenceHandlerInstance.decrementSequenceCount();
      if (numDisplay) {
        numDisplay.textContent = sequenceHandlerInstance.currentCount;
      }
    });
  }
});
