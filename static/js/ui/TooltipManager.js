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