// static/js/dom.js

/**
 * Manages sequence-related functionality
 */
console.log("✅ dom.js loaded! APP_STATE:", window.APP_STATE);

class SequenceHandler {
  constructor() {
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
      this.updateSequenceInputs(APP_STATE.testSeq || "1", true);
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
    console.log("🔄 Sending request to Flask:");
    console.log(" - numSequences:", numSequences);
    console.log(" - testingMode:", isTesting);
    console.log(" - testSeq:", isTesting ? APP_STATE.testSeq : "");

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
    document.querySelectorAll("#fileUpload, #numSequences").forEach(input => input.value = input.id === "numSequences" ? "1" : "");
    htmx.trigger(document.getElementById("numSequences"), "change");
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
  htmx.ajax("GET", "/get_sequence_inputs", {
    target: "#sequence-inputs-container",
    swap: "innerHTML",
    values: { numSequences: 1, testingMode: APP_STATE.testingMode }
  }).then(() => {
    initializeTabs();

    new SequenceHandler();
    new FormHandler();
    new TooltipManager();
  });

  // Ensure character counters update after form fields are pre-filled
  setTimeout(() => {
    document.querySelectorAll(".dynamic-sequence-input, #templateSequence").forEach(input => {
      if (input.value) {
        updateCharCount(input);
      }
    });
  }, 50);
  
  if (APP_STATE.testingMode) {
    document.getElementById("verbose_mode").checked = true;
  }
});
