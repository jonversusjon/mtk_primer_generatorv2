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