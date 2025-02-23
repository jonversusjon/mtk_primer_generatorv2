/**
* Handles horizontal scrolling for sequence tabs.
*/
export default class SequenceTabsScroller {
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