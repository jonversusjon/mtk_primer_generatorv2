import { useState, useEffect } from 'react';

/**
 * A hook for managing dark mode functionality
 * @returns {[boolean, Function]} A tuple containing current dark mode state and toggle function
 */
export const useDarkMode = () => {
  // Get stored preference or default to false (light mode)
  const getInitialMode = () => {
    const savedMode = localStorage.getItem('dark-mode');
    
    // Check for saved preference
    if (savedMode) {
      return savedMode === 'enabled';
    }
    
    // Check for system preference
    if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
      return true;
    }
    
    // Default to light mode
    return false;
  };
  
  const [darkMode, setDarkMode] = useState(getInitialMode);
  
  // Update localStorage and apply body class when darkMode changes
  useEffect(() => {
    // Store in localStorage
    localStorage.setItem('dark-mode', darkMode ? 'enabled' : 'disabled');
    
    // Apply to body class for global styling
    document.body.classList.toggle('dark-mode', darkMode);
  }, [darkMode]);
  
  // Toggle function
  const toggleDarkMode = () => {
    setDarkMode(!darkMode);
  };
  
  return [darkMode, toggleDarkMode];
};