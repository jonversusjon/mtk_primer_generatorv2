document.addEventListener("DOMContentLoaded", () => {
    // Parse data from the hidden JSON script tags
    const testSeqElement = document.getElementById("test-seq-data");
    const defaultTestSeq = testSeqElement ? JSON.parse(testSeqElement.textContent) : [];
    const testTemplateSeqElement = document.getElementById("test-template-seq");
    const defaultTestTemplateSeq = testTemplateSeqElement ? JSON.parse(testTemplateSeqElement.textContent) : "";
    const testingModeElement = document.getElementById("testing-mode-data");
    const TESTING_MODE = testingModeElement ? JSON.parse(testingModeElement.textContent) : false;
  
    // Pre-fill the template sequence if in testing mode
    const templateSequenceElement = document.getElementById("templateSequence");
    if (TESTING_MODE && templateSequenceElement && defaultTestTemplateSeq) {
      templateSequenceElement.value = defaultTestTemplateSeq;
    }
  
    /**
     * Updates the character count display for a given input or textarea.
     * @param {HTMLElement} el - The input or textarea element.
     */
    const updateCharCount = (el) => {
      const label = el.nextElementSibling;
      if (!label || !label.classList.contains("char-count-label")) return;
      const len = el.value.length;
      if (len > 0) {
        label.style.display = "block";
        label.textContent = `Length: ${len} bp`;
      } else {
        label.style.display = "none";
      }
    };
  
    /**
     * Attaches event listeners to all dynamic sequence inputs so that they update their char count.
     */
    const attachSequenceInputListeners = () => {
      document.querySelectorAll(".dynamic-sequence-input").forEach((input) => {
        input.addEventListener("input", () => updateCharCount(input));
        if (input.value) updateCharCount(input);
      });
    };
  
    /**
     * Fetches the HTML for the sequence inputs from the server via a Jinja-rendered partial.
     * Assumes your Flask route for `/get_sequence_inputs` renders a partial template (e.g. `sequence_input_tabs.html`)
     * that includes the complete structure (tabs, input groups, etc.).
     * @param {number|string} num - The number of sequence inputs to display.
     */
    const updateSequenceInputs = (num) => {
      const container = document.getElementById("sequence-container");
      // Clamp the number between 1 and 10
      num = Math.max(1, Math.min(10, Number(num)));
      fetch(`/get_sequence_inputs?numSequences=${num}`)
        .then((response) => response.text())
        .then((html) => {
          container.innerHTML = html;
          // After inserting the server-rendered HTML, attach listeners to the new inputs.
          attachSequenceInputListeners();
        })
        .catch((err) => console.error("Error fetching sequence inputs:", err));
    };
  
    /**
     * Clears the form by resetting the file input, sequence count (and thus the sequence inputs),
     * and any results.
     */
    const clearForm = () => {
      const fileUploadEl = document.getElementById("fileUpload");
      if (fileUploadEl) fileUploadEl.value = "";
      const numSequencesInput = document.getElementById("numSequences");
      if (numSequencesInput) {
        numSequencesInput.value = 1;
        updateSequenceInputs(1);
      }
      const resultsEl = document.getElementById("results");
      if (resultsEl) resultsEl.innerHTML = "";
    };
  
    /**
     * Switches the active tab in the Template Sequence section.
     * Designed to be called from inline HTML (via an onclick attribute).
     * @param {string} tabName - The ID of the tab pane to activate.
     * @param {Event} event - The event object from the click.
     */
    const openTab = (tabName, event) => {
      document.querySelectorAll(".tab-pane").forEach((pane) => pane.classList.remove("active"));
      document.querySelectorAll(".tab-button").forEach((button) => button.classList.remove("active"));
      const targetPane = document.getElementById(tabName);
      if (targetPane) targetPane.classList.add("active");
      if (event && event.currentTarget) event.currentTarget.classList.add("active");
    };
  
    /**
     * Switches the active sequence tab and its corresponding content pane.
     * @param {number} index - The index of the tab to activate.
     */
    const openSequenceTab = (index) => {
      const tabs = document.querySelectorAll(".sequence-tab");
      const panes = document.querySelectorAll(".sequence-pane");
      tabs.forEach(tab => tab.classList.remove("active"));
      panes.forEach(pane => pane.classList.remove("active"));
      if (tabs[index]) tabs[index].classList.add("active");
      if (panes[index]) panes[index].classList.add("active");
    };
  
    // Expose clearForm, openTab, and openSequenceTab globally so that inline HTML event handlers work
    window.clearForm = clearForm;
    window.openTab = openTab;
    window.openSequenceTab = openSequenceTab;
  
    // Initialize the sequence inputs on page load
    const numSequencesInput = document.getElementById("numSequences");
    if (numSequencesInput) {
      updateSequenceInputs(numSequencesInput.value);
      numSequencesInput.addEventListener("input", function () {
        updateSequenceInputs(this.value);
      });
    }
  
    // For the template sequence element, update char count on load and on input changes
    if (templateSequenceElement) {
      if (templateSequenceElement.value) updateCharCount(templateSequenceElement);
      templateSequenceElement.addEventListener("input", () => updateCharCount(templateSequenceElement));
    }
  });
  