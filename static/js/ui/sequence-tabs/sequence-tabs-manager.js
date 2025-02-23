import SequenceTabsScroller from "./sequence-tabs-scroller.js";

/**
 * Manages sequence-related functionality, including input handling, UI updates, and tab navigation.
 */
export default class SequenceTabsManager {
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
    * Initializes tab click behavior.
    */
    initializeTabs() {
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
     * Initializes tab scrolling for better navigation.
     */
    initializeTabScroller(containerId) {
        new SequenceTabsScroller(containerId);
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

    /**
 * Increments the sequence count up to a maximum of 10.
 * - Makes the next tab’s nav button visible.
 * - Activates the newly added tab.
 */
    incrementSequenceCount() {
        if (this.currentCount < 10) {
            this.currentCount++;
            // Update the nav buttons and tab contents so that tabs with index > currentCount are hidden.
            this.updateSequenceInputs(this.currentCount);
            // Activate the newly added tab (only one tab active at a time)
            const newTabBtn = document.querySelector(`.sequence-tab-btn[data-tab-index="${this.currentCount}"]`);
            if (newTabBtn) {
                newTabBtn.click();
            }
        }
    }

    /**
     * Decrements the sequence count down to a minimum of 1.
     * - Clears all input data on the current (nth) tab.
     * - Hides the current tab’s nav button.
     * - Activates the previous tab.
     */
    decrementSequenceCount() {
        if (this.currentCount > 1) {
            // Clear data on the current (nth) tab that is about to be removed.
            const tabContent = document.querySelector(`.sequence-tab-content[data-tab-index="${this.currentCount}"]`);
            if (tabContent) {
                // Clear the DNA sequence textarea.
                const seqTextarea = tabContent.querySelector("textarea.dynamic-sequence-input");
                if (seqTextarea) {
                    seqTextarea.value = "";
                    this.updateCharCount(seqTextarea);
                }
                // Clear the primer name input.
                const primerInput = tabContent.querySelector("input.primer-name-input");
                if (primerInput) {
                    primerInput.value = "";
                }
                // Reset the MTK part number select to its default (first option).
                const mtkSelect = tabContent.querySelector("select.mtk-part-select");
                if (mtkSelect) {
                    mtkSelect.selectedIndex = 0;
                }
            }
            this.currentCount--;
            // Update the UI to hide the removed tab.
            this.updateSequenceInputs(this.currentCount);
            // Activate the now last visible tab.
            const lastTabBtn = document.querySelector(`.sequence-tab-btn[data-tab-index="${this.currentCount}"]`);
            if (lastTabBtn) {
                lastTabBtn.click();
            }
        }
    }

}

// Initialize the sequence handler when the DOM is loaded
document.addEventListener("DOMContentLoaded", () => new SequenceTabsManager());
