/**
 * Manages sequence-related functionality, including input handling, UI updates, and tab navigation.
 */
class SequenceTabsHandler {
  constructor() {
      const numSequencesInput = document.getElementById("numSequencesInput");
      this.currentCount = numSequencesInput ? parseInt(numSequencesInput.value, 10) || 1 : 1;

      this.initializeState();
      this.initializeTabs();
      this.initializeTabScroller("sequence-tabs-container");
  }

  /**
   * Initializes the application state based on whether testing mode is enabled.
   */
  initializeState() {
      if (typeof APP_STATE !== "undefined" && APP_STATE.testingMode) {
          this.setTestDefaults();
          const numTestSeq = APP_STATE.testSeq?.length || 1;
          console.log("✅ Test mode detected. Setting sequence count to", numTestSeq);
          this.updateSequenceInputs(numTestSeq, true);
      } else {
          this.updateSequenceInputs(1, false);
      }
  }

  /**
   * Initializes tab scrolling for better navigation.
   */
  initializeTabScroller(containerId) {
      new TabScroller(containerId);
  }

  /**
   * Populates form fields with test data if test mode is enabled.
   */
  setTestDefaults() {
      const templateSequence = document.getElementById("templateSequence");
      if (templateSequence && APP_STATE.testTemplateSeq) {
          templateSequence.value = APP_STATE.testTemplateSeq;
          this.updateCharCount(templateSequence);
      }

      if (Array.isArray(APP_STATE.testSeq)) {
          APP_STATE.testSeq.forEach((testObj, index) => {
              const seqIndex = index + 1;
              const { sequence = "", primerName = `Primer ${seqIndex}`, mtkPartNum = "" } = testObj || {};

              const seqTextarea = document.getElementById(`sequenceInput${seqIndex}`);
              const primerInput = document.getElementById(`primerName${seqIndex}`);
              const mtkSelect = document.getElementById(`mtkPartNum${seqIndex}`);

              if (seqTextarea) {
                  seqTextarea.value = sequence;
                  this.updateCharCount(seqTextarea);
              }
              if (primerInput) primerInput.value = primerName;
              if (mtkSelect) mtkSelect.value = mtkPartNum || mtkSelect.selectedIndex;
          });
      }
  }

  /**
   * Updates the number of visible sequence tabs and synchronizes the UI.
   */
  updateVisibleTabs(newCount) {
      document.querySelectorAll(".sequence-tab-btn, .sequence-tab-content").forEach(el => {
          const idx = parseInt(el.getAttribute("data-tab-index"), 10);
          el.classList.toggle("hidden", idx > newCount);
      });

      const numInput = document.getElementById("numSequencesInput");
      if (numInput) numInput.value = newCount;
  }

  /**
   * Updates the displayed sequences and their corresponding UI elements.
   */
  updateSequenceInputs(newCount) {
      this.currentCount = parseInt(newCount, 10);
      this.updateVisibleTabs(this.currentCount);

      // Ensure at least one active tab
      const activeBtn = document.querySelector(".sequence-tab-btn.active");
      if (!activeBtn) {
          document.querySelector(".sequence-tab-btn:not(.hidden)")?.click();
      }
  }

  /**
   * Updates the character count label for a given input.
   */
  updateCharCount(element) {
      const label = element.nextElementSibling;
      if (label?.classList.contains("char-count-label")) {
          label.textContent = `Length: ${element.value.length} bp`;
          label.style.display = element.value.length ? "block" : "none";
      }
  }
}

/**
* Initializes tab click behavior.
*/
export function initializeTabs() {
  document.querySelectorAll(".sequence-tab-btn").forEach(btn => {
      btn.addEventListener("click", () => {
          const tabIndex = btn.getAttribute("data-tab-index");

          document.querySelectorAll(".sequence-tab-btn").forEach(el => el.classList.remove("active"));
          btn.classList.add("active");

          document.querySelectorAll(".sequence-tab-content").forEach(el => {
              el.classList.toggle("active", el.getAttribute("data-tab-index") === tabIndex);
              el.style.opacity = el.classList.contains("active") ? "1" : "0";
          });
      });
  });
}

/**
* Handles horizontal scrolling for sequence tabs.
*/
class TabScroller {
  constructor(containerId) {
      this.container = document.getElementById(containerId);
      if (!this.container) return;

      this.tabList = this.container.querySelector(".sequence-tabs-nav");
      this.setupScrollButtons();
      this.checkScrollButtons();

      window.addEventListener("resize", () => this.checkScrollButtons());
      this.tabList.addEventListener("scroll", () => this.checkScrollButtons());
  }

  setupScrollButtons() {
      this.leftButton = this.createButton("‹", "left");
      this.rightButton = this.createButton("›", "right");

      this.container.insertBefore(this.leftButton, this.tabList);
      this.container.appendChild(this.rightButton);
  }

  createButton(text, direction) {
      const button = document.createElement("button");
      button.className = `tab-scroll-button ${direction}`;
      button.textContent = text;
      button.addEventListener("click", () => this.scroll(direction));
      return button;
  }

  checkScrollButtons() {
      const { scrollLeft, scrollWidth, clientWidth } = this.tabList;
      this.leftButton.style.display = scrollLeft > 0 ? "flex" : "none";
      this.rightButton.style.display = scrollLeft < (scrollWidth - clientWidth - 1) ? "flex" : "none";
  }

  scroll(direction) {
      this.tabList.scrollBy({
          left: (direction === "left" ? -1 : 1) * this.tabList.clientWidth * 0.5,
          behavior: "smooth"
      });
  }
}

// Initialize the sequence handler when the DOM is loaded
document.addEventListener("DOMContentLoaded", () => new SequenceTabsHandler());
