// static/js/dom.js

/**
 * Manages sequence-related functionality
 */
console.log("✅ dom.js loaded! APP_STATE:", window.APP_STATE);

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

function initializeTabs() {
  document.querySelectorAll(".sequence-tab-btn").forEach(btn => {
    btn.addEventListener("click", () => {
      const tabIndex = btn.getAttribute("data-tab-index");

      // Deactivate all nav buttons, then activate the clicked one.
      document.querySelectorAll(".sequence-tab-btn").forEach(el => el.classList.remove("active"));
      btn.classList.add("active");

      // Hide all tab contents.
      document.querySelectorAll(".sequence-tab-content").forEach(el => {
        el.classList.remove("active");
        el.style.opacity = "0";
      });

      // Activate the content pane matching the clicked nav button.
      const activeContent = document.querySelector(`.sequence-tab-content[data-tab-index="${tabIndex}"]`);
      if (activeContent) {
        activeContent.classList.add("active");
        setTimeout(() => {
          activeContent.style.opacity = "1";
        }, 50);
      }
    });
  });

  // Initialize TabScroller for horizontal scrolling, if needed.
  new TabScroller("sequence-tabs-container");
}

class TabScroller {
  constructor(containerId) {
    this.container = document.getElementById(containerId);
    if (!this.container) return;

    this.tabList = this.container.querySelector('.sequence-tabs-nav');
    this.setupScrollButtons();
    this.checkScrollButtons();

    window.addEventListener('resize', () => this.checkScrollButtons());
    this.tabList.addEventListener('scroll', () => this.checkScrollButtons());
  }

  setupScrollButtons() {
    const leftButton = document.createElement('button');
    leftButton.className = 'tab-scroll-button left';
    leftButton.innerHTML = '‹';
    leftButton.addEventListener('click', () => this.scroll('left'));

    const rightButton = document.createElement('button');
    rightButton.className = 'tab-scroll-button right';
    rightButton.innerHTML = '›';
    rightButton.addEventListener('click', () => this.scroll('right'));

    this.container.insertBefore(leftButton, this.tabList);
    this.container.appendChild(rightButton);

    this.leftButton = leftButton;
    this.rightButton = rightButton;
  }

  checkScrollButtons() {
    const { scrollLeft, scrollWidth, clientWidth } = this.tabList;
    this.leftButton.style.display = scrollLeft > 0 ? 'flex' : 'none';
    this.rightButton.style.display = scrollLeft < (scrollWidth - clientWidth - 1) ? 'flex' : 'none';
  }

  scroll(direction) {
    const scrollAmount = this.tabList.clientWidth * 0.5;
    const newScrollPosition = this.tabList.scrollLeft +
      (direction === 'left' ? -scrollAmount : scrollAmount);

    this.tabList.scrollTo({
      left: newScrollPosition,
      behavior: 'smooth'
    });
  }
}

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
