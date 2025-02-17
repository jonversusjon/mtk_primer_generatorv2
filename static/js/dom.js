// Initialize state and constants from embedded data
const APP_STATE = {
  testSeq: JSON.parse(document.getElementById("test-seq-data")?.textContent || "[]"),
  testTemplateSeq: JSON.parse(document.getElementById("test-template-seq")?.textContent || "\"\""),
  testingMode: JSON.parse(document.getElementById("testing-mode-data")?.textContent || "false"),
};

class SequenceManager {
  constructor() {
    this.container = document.getElementById("sequence-container");
    this.numSequencesInput = document.getElementById("numSequences");
    this.templateSequence = document.getElementById("templateSequence");
    this.setupEventListeners();
    this.initializeState();
  }

  setupEventListeners() {
    if (this.numSequencesInput) {
      this.numSequencesInput.addEventListener("input", (e) => this.updateSequenceInputs(e.target.value));
    }
    if (this.templateSequence) {
      this.templateSequence.addEventListener("input", () => this.updateCharCount(this.templateSequence));
    }
  }

  initializeState() {
    // Initialize with current number of sequences
    if (this.numSequencesInput) {
      this.updateSequenceInputs(this.numSequencesInput.value);
    }
    // Initialize template sequence if in testing mode
    if (APP_STATE.testingMode && this.templateSequence && APP_STATE.testTemplateSeq) {
      this.templateSequence.value = APP_STATE.testTemplateSeq;
      this.updateCharCount(this.templateSequence);
    }
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

  attachSequenceInputListeners() {
    document.querySelectorAll(".dynamic-sequence-input").forEach(input => {
      input.addEventListener("input", () => this.updateCharCount(input));
      if (input.value) this.updateCharCount(input);
    });
  }

  async updateSequenceInputs(num) {
    num = Math.max(1, Math.min(10, Number(num)));
    try {
      const response = await fetch(`/get_sequence_inputs?numSequences=${num}`);
      const html = await response.text();
      if (this.container) {
        this.container.innerHTML = html;
        this.attachSequenceInputListeners();
        this.initializeTabs();
      }
    } catch (error) {
      console.error("Error updating sequence inputs:", error);
    }
  }

  initializeTabs() {
    const tabButtons = document.querySelectorAll('.sequence-tab-btn');
    tabButtons.forEach(button => {
      button.addEventListener('click', () => this.switchTab(button.dataset.tabIndex));
    });
    // Activate first tab
    if (tabButtons.length > 0) {
      this.switchTab(0);
    }
  }

  switchTab(index) {
    document.querySelectorAll('.sequence-tab-btn').forEach(btn => {
      btn.classList.toggle('active', btn.dataset.tabIndex === index.toString());
    });
    document.querySelectorAll('.sequence-tab-content').forEach(content => {
      content.classList.toggle('active', content.dataset.tabIndex === index.toString());
    });
  }
}

class FormManager {
  constructor() {
    this.runDesignButton = document.getElementById("run-design-primer");
    this.clearFormButton = document.getElementById("clear-form");
    this.resultsElement = document.getElementById("results");
    this.setupEventListeners();
  }

  setupEventListeners() {
    if (this.clearFormButton) {
      this.clearFormButton.addEventListener("click", () => this.clearForm());
    }
    if (this.runDesignButton) {
      this.runDesignButton.addEventListener("click", () => this.handleDesignSubmission());
    }
  }

  clearForm() {
    document.getElementById("fileUpload")?.value = "";
    const numSequencesInput = document.getElementById("numSequences");
    if (numSequencesInput) {
      numSequencesInput.value = "1";
      numSequencesInput.dispatchEvent(new Event("input"));
    }
    if (this.resultsElement) {
      this.resultsElement.innerHTML = "";
    }
  }

  async handleDesignSubmission() {
    if (!this.runDesignButton || !this.resultsElement) return;

    this.runDesignButton.disabled = true;
    this.clearFormButton.disabled = true;
    this.resultsElement.innerHTML = "<p>Designing primers...</p>";

    try {
      // Replace with actual submission logic
      await new Promise(resolve => setTimeout(resolve, 3000));
      this.resultsElement.innerHTML = "<p>Primer design completed.</p>";
    } catch (error) {
      this.resultsElement.innerHTML = "<p>Error designing primers.</p>";
      console.error("Error:", error);
    } finally {
      this.runDesignButton.disabled = false;
      this.clearFormButton.disabled = false;
    }
  }
}

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

// Initialize everything when DOM is ready
document.addEventListener("DOMContentLoaded", () => {
  const sequenceManager = new SequenceManager();
  const formManager = new FormManager();
  const tooltipManager = new TooltipManager();

  // Set testing mode checkbox if needed
  if (APP_STATE.testingMode) {
    document.getElementById("verbose_mode").checked = true;
  }
});