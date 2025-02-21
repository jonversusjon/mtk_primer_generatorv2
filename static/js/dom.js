// static/js/dom.js

/**
 * Manages sequence-related functionality
 */
console.log("✅ dom.js loaded! APP_STATE:", window.APP_STATE);

class SequenceHandler {
  constructor() {
    this.currentCount = 1; /* this is the input_sequences counter */
    this.setupCharacterCounters();
    this.initializeState();
  }

  setupCharacterCounters() {
    document.querySelectorAll(".dynamic-sequence-input, #templateSequence").forEach(input => {
      input.addEventListener("input", () => this.updateCharCount(input));
      if (input.value) this.updateCharCount(input);
    });
  }

  initializeState() {
    if (APP_STATE.testingMode) {
      this.setTestDefaults();
      const numTestSeq = APP_STATE.testSeq.length > 0 ? APP_STATE.testSeq.length : 1;
      console.log("✅ Test mode detected. Setting sequence count to", numTestSeq);

      this.updateSequenceInputs(numTestSeq, true);
    } else {
      this.updateSequenceInputs("1", false);
    }
  }

  setTestDefaults() {
    const templateSequence = document.getElementById("templateSequence");
    if (templateSequence) {
      templateSequence.value = APP_STATE.testTemplateSeq;
      this.updateCharCount(templateSequence);
    }
  }


  updateSequenceInputs(numSequences, isTesting) {

    htmx.ajax("GET", "/get_sequence_inputs", {
      target: "#sequence-inputs-container",
      swap: "innerHTML",
      values: {
        numSequences: numSequences,
        testingMode: isTesting,
        testSeq: isTesting ? APP_STATE.testSeq : "",
      }
    }).then(() => {
      console.log("✅ Flask response received, updating UI");
      initializeTabs();
      this.setupCharacterCounters();
    });
  }

  updateCharCount(element) {
    const label = element.nextElementSibling;
    if (!label?.classList.contains("char-count-label")) return;

    label.style.display = element.value.length ? "block" : "none";
    label.textContent = `Length: ${element.value.length} bp`;
  }

  changeSequenceCount(delta) {
    const newCount = this.currentCount + delta;
    console.log("Attempting to change sequence count from", this.currentCount, "to", newCount);
    if (newCount < 1 || newCount > 10) {
      console.log("New count out of bounds, ignoring.");
      return;
    }
    this.currentCount = newCount;
    this.updateSequenceInputs(this.currentCount, APP_STATE.testingMode);
  }

}

/**
 * Manages form clearing and reset functionality
 */
class FormHandler {
  constructor() {
    this.clearFormButton = document.getElementById("clear-form");
    if (this.clearFormButton) {
      this.clearFormButton.addEventListener("click", () => this.clearForm());
    }
  }

  clearForm() {
    // Clear file upload fields
    document.querySelectorAll("#fileUpload").forEach(input => input.value = "");
    // Reset the sequence count (if our SequenceHandler instance exists)
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

/**
 * Manages tooltip functionality
 */
class TooltipManager {
  constructor() {
    this.setupTooltips();
  }

  setupTooltips() {
    document.querySelectorAll(".info-icon-container").forEach(container => {
      const tooltip = container.querySelector(".info-tooltip");
      const icon = container.querySelector(".info-icon");

      if (icon && tooltip) {
        icon.addEventListener("click", () => this.toggleTooltip(tooltip));
        container.addEventListener("mouseenter", () => this.positionTooltip(tooltip));
      }
    });

    document.addEventListener("click", event => {
      if (!event.target.closest(".info-icon-container")) {
        document.querySelectorAll(".info-tooltip").forEach(tooltip => tooltip.style.display = "none");
      }
    });
  }

  toggleTooltip(tooltip) {
    tooltip.style.display = tooltip.style.display === "block" ? "none" : "block";
    if (tooltip.style.display === "block") this.positionTooltip(tooltip);
  }

  positionTooltip(tooltip) {
    const rect = tooltip.getBoundingClientRect();
    tooltip.style.left = rect.right > window.innerWidth ? "auto" : "50%";
    tooltip.style.right = rect.right > window.innerWidth ? "0" : "auto";
    tooltip.style.transform = rect.right > window.innerWidth ? "none" : "translateX(-50%)";
    tooltip.style.top = rect.bottom > window.innerHeight ? "125%" : "auto";
    tooltip.style.bottom = rect.bottom > window.innerHeight ? "auto" : "125%";
  }
}

/**
 * Initializes tab functionality for sequence inputs
 */
function initializeTabs() {
  document.querySelectorAll(".sequence-tab-btn").forEach(btn => {
    btn.addEventListener("click", () => {
      const tabIndex = btn.getAttribute("data-tab-index");

      document.querySelectorAll(".sequence-tab-btn, .sequence-tab-content").forEach(el => el.classList.remove("active"));
      btn.classList.add("active");
      document.querySelector(`.sequence-tab-content[data-tab-index="${tabIndex}"]`)?.classList.add("active");
    });
  });
}

// Initialize everything when DOM is ready
document.addEventListener("DOMContentLoaded", () => {
  // Instantiate the SequenceHandler and store it globally for access (e.g., from the clearForm method or button events)
  window.sequenceHandlerInstance = new SequenceHandler();
  new FormHandler();
  new TooltipManager();

  // Ensure character counters update after form fields are pre-filled
  setTimeout(() => {
    document.querySelectorAll(".dynamic-sequence-input, #templateSequence").forEach(input => {
      if (input.value) {
        window.sequenceHandlerInstance.updateCharCount(input);
      }
    });
  }, 50);

  if (APP_STATE.testingMode) {
    document.getElementById("verbose_mode").checked = true;
  }

  const numInput = document.getElementById("numSequencesInput");
  const numDisplay = document.getElementById("numDisplay");
  const incrementBtn = document.getElementById("incrementBtn");
  const decrementBtn = document.getElementById("decrementBtn");

  if (numInput && numDisplay && incrementBtn && decrementBtn) {
    function updateSequenceCount(delta) {
      let currentValue = parseInt(numInput.value, 10) || 0;
      let newValue = currentValue + delta;
      if (newValue < 1) newValue = 1;
      if (newValue > 10) newValue = 10;
      numInput.value = newValue;
      numDisplay.textContent = newValue;

      if (window.sequenceHandlerInstance) {
        window.sequenceHandlerInstance.updateSequenceInputs(newValue, APP_STATE.testingMode);
      } else {
        console.log("can't update number of sequence inputs without a sequenceHandlerInstance.")
      }
    }

    incrementBtn.addEventListener("click", () => updateSequenceCount(1));
    decrementBtn.addEventListener("click", () => updateSequenceCount(-1));
  } else {
    console.warn("Sequence input or buttons not found in the DOM.");
  }
});
