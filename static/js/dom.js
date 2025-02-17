// static/js/dom.js


/**
 * Manages sequence-related functionality
 */
class SequenceHandler {
  constructor() {
    this.setupCharacterCounters();
    this.initializeState();
  }

  setupCharacterCounters() {
    // Setup for template sequence
    const templateSequence = document.getElementById("templateSequence");
    if (templateSequence) {
      templateSequence.addEventListener("input", () => this.updateCharCount(templateSequence));
      if (templateSequence.value) this.updateCharCount(templateSequence);
    }

    // Setup for dynamic sequence inputs
    document.querySelectorAll(".dynamic-sequence-input").forEach(input => {
      input.addEventListener("input", () => this.updateCharCount(input));
      if (input.value) this.updateCharCount(input);
    });
  }

  initializeState() {
    if (APP_STATE.testingMode) {
      // Initialize template sequence if it exists
      const templateSequence = document.getElementById("templateSequence");
      if (templateSequence && APP_STATE.testTemplateSeq) {
        templateSequence.value = APP_STATE.testTemplateSeq;
        this.updateCharCount(templateSequence);
      }

      // Update sequence inputs with test data
      this.updateSequenceInputsWithTestData();
    } else {
      // If not in testing mode, load without test data
      this.updateSequenceInputs("1");
    }
  }

  updateSequenceInputsWithTestData() {
    const numSequences = document.getElementById("numSequences")?.value || "1";
    htmx.ajax('GET', '/get_sequence_inputs', {
      target: '#sequence-container',
      swap: 'innerHTML',
      values: {
        numSequences: numSequences,
        testingMode: true,
        testSeq: APP_STATE.testSeq
      }
    }).then(() => {
      initializeTabs();
      this.setupCharacterCounters();
    });
  }

  updateSequenceInputs(value) {
    const numSequences = parseInt(value) || 1;
    htmx.ajax('GET', '/get_sequence_inputs', {
      target: '#sequence-container',
      swap: 'innerHTML',
      values: {
        numSequences: numSequences,
        testingMode: false
      }
    }).then(() => {
      initializeTabs();
      this.setupCharacterCounters();
    });
  }

  updateCharCount(element) {
    const label = element.nextElementSibling;
    if (!label?.classList.contains("char-count-label")) return;

    const length = element.value.length;
    label.style.display = length > 0 ? "block" : "none";
    if (length > 0) {
      label.textContent = `Length: ${length} bp`;
    }
  }
}

/**
 * Manages form clearing and reset functionality
 */
class FormHandler {
  constructor() {
    this.clearFormButton = document.getElementById("clear-form");
    this.setupEventListeners();
  }

  setupEventListeners() {
    if (this.clearFormButton) {
      this.clearFormButton.addEventListener("click", () => this.clearForm());
    }
  }

  clearForm() {
    // Reset file upload
    const fileUpload = document.getElementById("fileUpload");
    if (fileUpload) {
        fileUpload.value = "";
    }

    // Reset number of sequences
    const numSequencesInput = document.getElementById("numSequences");
    if (numSequencesInput) {
      numSequencesInput.value = "1";
      // Trigger HTMX request to update sequence inputs
      htmx.trigger(numSequencesInput, 'change');
    }

    // Clear results
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
    this.setupGlobalListeners();
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
  }

  setupGlobalListeners() {
    document.addEventListener("click", (event) => {
      if (!event.target.closest(".info-icon-container")) {
        document.querySelectorAll(".info-tooltip").forEach(tooltip => {
          tooltip.style.display = "none";
        });
      }
    });
  }

  toggleTooltip(tooltip) {
    tooltip.style.display = tooltip.style.display === "block" ? "none" : "block";
    if (tooltip.style.display === "block") {
      this.positionTooltip(tooltip);
    }
  }

  positionTooltip(tooltip) {
    const rect = tooltip.getBoundingClientRect();
    const isOffScreenRight = rect.right > window.innerWidth;
    const isOffScreenBottom = rect.bottom > window.innerHeight;

    tooltip.style.left = isOffScreenRight ? "auto" : "50%";
    tooltip.style.right = isOffScreenRight ? "0" : "auto";
    tooltip.style.transform = isOffScreenRight ? "none" : "translateX(-50%)";
    tooltip.style.top = isOffScreenBottom ? "125%" : "auto";
    tooltip.style.bottom = isOffScreenBottom ? "auto" : "125%";
  }
}

/**
 * Updates sequence inputs based on the selected number of sequences
 * @param {string|number} value - The number of sequences to generate
 */
function updateSequenceInputs(value) {
  const numSequences = parseInt(value) || 1;
  const isTestingMode = APP_STATE.testingMode;

  htmx.ajax('GET', '/get_sequence_inputs', {
    target: '#sequence-container',
    swap: 'innerHTML',
    values: {
      numSequences: numSequences,
      testingMode: isTestingMode
    }
  }).then(() => {
    initializeTabs();
    const sequenceHandler = new SequenceHandler();
  });
}

/**
 * Initializes tab functionality for sequence inputs
 */
function initializeTabs() {
  const tabBtns = document.querySelectorAll('.sequence-tab-btn');
  const tabContents = document.querySelectorAll('.sequence-tab-content');

  tabBtns.forEach(btn => {
    btn.addEventListener('click', () => {
      const tabIndex = btn.getAttribute('data-tab-index');

      // Remove active class from all buttons and contents
      tabBtns.forEach(b => b.classList.remove('active'));
      tabContents.forEach(c => c.classList.remove('active'));

      // Add active class to clicked button and corresponding content
      btn.classList.add('active');
      document.querySelector(`.sequence-tab-content[data-tab-index="${tabIndex}"]`)
        ?.classList.add('active');
    });
  });
}

// Initialize everything when DOM is ready
document.addEventListener("DOMContentLoaded", () => {
  console.log("DOM Content Loaded. Testing mode:", APP_STATE.testingMode);

  // Initial load of sequence input
  htmx.ajax('GET', '/get_sequence_inputs', {
    target: '#sequence-container',
    swap: 'innerHTML',
    values: {
      numSequences: 1,
      testingMode: APP_STATE.testingMode  // Add this line
    }
  }).then(() => {
    initializeTabs();
    const sequenceHandler = new SequenceHandler();
    const formHandler = new FormHandler();
    const tooltipManager = new TooltipManager();
  });

  // Set testing mode checkbox if needed
  if (APP_STATE.testingMode) {
    document.getElementById("verbose_mode").checked = true;
  }
});