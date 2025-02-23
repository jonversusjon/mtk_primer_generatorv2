/**
 * Manages tooltip functionality
 */
class TooltipManager {
  constructor() {
      this.setupTooltips();
  }

  /**
   * Attaches event listeners to tooltip elements.
   */
  setupTooltips() {
      document.querySelectorAll(".info-icon-container").forEach(container => {
          const tooltip = container.querySelector(".info-tooltip");
          const icon = container.querySelector(".info-icon");

          if (icon && tooltip) {
              icon.addEventListener("click", (event) => {
                  event.stopPropagation();
                  this.toggleTooltip(tooltip);
              });

              container.addEventListener("mouseenter", () => this.showTooltip(tooltip));
              container.addEventListener("mouseleave", () => this.hideTooltip(tooltip));
          }
      });

      // Hide tooltips if user clicks anywhere outside a tooltip
      document.addEventListener("click", () => {
          this.hideAllTooltips();
      });
  }

  /**
   * Shows the tooltip and positions it correctly.
   */
  showTooltip(tooltip) {
      tooltip.style.visibility = "hidden";
      tooltip.style.display = "block";
      this.positionTooltip(tooltip);
      tooltip.style.visibility = "visible";
      tooltip.setAttribute("aria-hidden", "false");
  }

  /**
   * Hides the tooltip.
   */
  hideTooltip(tooltip) {
      tooltip.style.display = "none";
      tooltip.setAttribute("aria-hidden", "true");
  }

  /**
   * Hides all tooltips.
   */
  hideAllTooltips() {
      document.querySelectorAll(".info-tooltip").forEach(tooltip => this.hideTooltip(tooltip));
  }

  /**
   * Toggles the visibility of a tooltip.
   */
  toggleTooltip(tooltip) {
      if (tooltip.style.display === "block") {
          this.hideTooltip(tooltip);
      } else {
          this.showTooltip(tooltip);
      }
  }

  /**
   * Positions the tooltip to ensure it does not overflow the viewport.
   */
  positionTooltip(tooltip) {
      const icon = tooltip.closest(".info-icon-container")?.querySelector(".info-icon");
      if (!icon) return;

      const tooltipRect = tooltip.getBoundingClientRect();
      const iconRect = icon.getBoundingClientRect();

      // Default positioning
      let top = iconRect.bottom + window.scrollY + 6; // Add spacing below the icon
      let left = iconRect.left + window.scrollX + iconRect.width / 2 - tooltipRect.width / 2;

      // Prevent tooltip from going off-screen horizontally
      if (left < 10) left = 10; // Prevent left overflow
      if (left + tooltipRect.width > window.innerWidth - 10) {
          left = window.innerWidth - tooltipRect.width - 10; // Prevent right overflow
      }

      // Prevent tooltip from going off-screen vertically
      if (tooltipRect.bottom > window.innerHeight) {
          top = iconRect.top + window.scrollY - tooltipRect.height - 6; // Position above the icon
      }

      // Apply styles
      tooltip.style.left = `${left}px`;
      tooltip.style.top = `${top}px`;
  }
}

// Initialize when DOM is loaded
document.addEventListener("DOMContentLoaded", () => {
  new TooltipManager();
});