/**
 * Manages sequence-related functionality
 */

class SequenceHandler {
    constructor() {
      const initialCount = parseInt(document.getElementById("numSequencesInput").value, 10) || 1;
      this.currentCount = initialCount;
      this.setupCharacterCounters();
      this.initializeState();
    }
  
    setupCharacterCounters() {
      // Use single class for all sequence textareas
      document.querySelectorAll(".dynamic-sequence-input").forEach(input => {
        input.addEventListener("input", () => this.updateCharCount(input));
        if (input.value) this.updateCharCount(input);
      });
    }
  
    initializeState() {
      if (APP_STATE.testingMode) {
        this.setTestDefaults();
        const numTestSeq = (APP_STATE.testSeq && APP_STATE.testSeq.length) ? APP_STATE.testSeq.length : 1;
        console.log("✅ Test mode detected. Setting sequence count to", numTestSeq);
        this.updateSequenceInputs(numTestSeq, true);
      } else {
        this.updateSequenceInputs(1, false);
      }
    }
  
    setTestDefaults() {
      const templateSequence = document.getElementById("templateSequence");
      if (templateSequence) {
        templateSequence.value = APP_STATE.testTemplateSeq;
        this.updateCharCount(templateSequence);
      }
  
      if (APP_STATE.testSeq && Array.isArray(APP_STATE.testSeq)) {
        APP_STATE.testSeq.forEach((testObj, index) => {
          const seqIndex = index + 1;
          const seqTextarea = document.getElementById(`sequenceInput${seqIndex}`);
          const primerInput = document.getElementById(`primerName${seqIndex}`);
          const mtkSelect = document.getElementById(`mtkPartNum${seqIndex}`);
  
          let sequenceValue = "";
          if (typeof testObj === "string") {
            sequenceValue = testObj;
          } else if (typeof testObj === "object" && testObj.sequence) {
            sequenceValue = testObj.sequence;
          }
  
          if (seqTextarea) {
            seqTextarea.value = sequenceValue;
            this.updateCharCount(seqTextarea);
          }
  
          if (primerInput) {
            if (typeof testObj === "object" && testObj.primerName) {
              primerInput.value = testObj.primerName;
            } else {
              primerInput.value = `Primer ${seqIndex}`;
            }
          }
  
          if (mtkSelect) {
            if (typeof testObj === "object" && testObj.mtkPartNum) {
              mtkSelect.value = testObj.mtkPartNum;
            } else {
              mtkSelect.selectedIndex = 0;
            }
          }
        });
      }
    }
  
    updateVisibleTabs(newCount) {
      document.querySelectorAll(".sequence-tab-btn").forEach(btn => {
        const idx = parseInt(btn.getAttribute("data-tab-index"), 10);
        if (idx <= newCount) {
          btn.classList.remove("hidden");
        } else {
          btn.classList.add("hidden");
        }
      });
  
      document.querySelectorAll(".sequence-tab-content").forEach(content => {
        const idx = parseInt(content.getAttribute("data-tab-index"), 10);
        if (idx <= newCount) {
          content.classList.remove("hidden");
        } else {
          content.classList.add("hidden");
        }
      });
  
      const numInput = document.getElementById("numSequencesInput");
      if (numInput) numInput.value = newCount;
    }
  
    updateSequenceInputs(newCount, isTesting) {
      this.currentCount = parseInt(newCount, 10);
      this.updateVisibleTabs(this.currentCount);
  
      const activeBtn = document.querySelector(".sequence-tab-btn.active");
      if (!activeBtn) {
        const firstVisibleBtn = document.querySelector(".sequence-tab-btn:not(.hidden)");
        if (firstVisibleBtn) {
          firstVisibleBtn.classList.add("active");
          firstVisibleBtn.click();
        }
      }
    }
  
    updateCharCount(element) {
      const label = element.nextElementSibling;
      if (!label || !label.classList.contains("char-count-label")) return;
      label.style.display = element.value.length ? "block" : "none";
      label.textContent = `Length: ${element.value.length} bp`;
    }
  
    incrementSequenceCount() {
      if (this.currentCount < 10) {
        this.currentCount++;
        this.updateVisibleTabs(this.currentCount);
        const newTabButton = document.querySelector(`.sequence-tab-btn[data-tab-index="${this.currentCount}"]`);
        if (newTabButton) {
          newTabButton.click();
        }
      }
    }
  
    decrementSequenceCount() {
      if (this.currentCount > 1) {
        this.clearTabData(this.currentCount);
        this.currentCount--;
        this.updateVisibleTabs(this.currentCount);
        const newTabButton = document.querySelector(`.sequence-tab-btn[data-tab-index="${this.currentCount}"]`);
        if (newTabButton) {
          newTabButton.click();
        }
      }
    }
  
    clearTabData(tabIndex) {
      const tabContainer = document.querySelector(`.sequence-tab-content[data-tab-index="${tabIndex}"]`);
      if (tabContainer) {
        tabContainer.querySelectorAll(".error-message").forEach(el => {
          el.innerHTML = "";
        });
      }
  
      const seqTextarea = document.getElementById(`sequenceInput${tabIndex}`);
      const primerInput = document.getElementById(`primerName${tabIndex}`);
      const mtkSelect = document.getElementById(`mtkPartNum${tabIndex}`);
  
      if (seqTextarea) {
        seqTextarea.value = "";
        this.updateCharCount(seqTextarea);
      }
      if (primerInput) {
        primerInput.value = "";
      }
      if (mtkSelect) {
        mtkSelect.selectedIndex = 0;
      }
    }
  }