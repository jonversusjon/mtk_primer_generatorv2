import { useState, useEffect } from "react";

function useDebouncedLocalStorage(key, initialValue, delay = 500) {
  // Initialize state from localStorage (or fallback to initialValue)
  const [value, setValue] = useState(() => {
    const storedValue = localStorage.getItem(key);
    return storedValue ? JSON.parse(storedValue) : initialValue;
  });

  useEffect(() => {
    // Set up a timeout to update localStorage after the delay
    const handler = setTimeout(() => {
      localStorage.setItem(key, JSON.stringify(value));
    }, delay);

    // Clear the timeout if value changes before delay completes
    return () => {
      clearTimeout(handler);
    };
  }, [value, key, delay]);

  return [value, setValue];
}

export default useDebouncedLocalStorage;
