import React from "react";
import { useMediaQuery } from "react-responsive";
import { FaExclamationCircle } from "react-icons/fa";

export default function ErrorBubble({ errorCount, onClick }) {
  const isMobile = useMediaQuery({ maxWidth: 859 });

  if (!isMobile || errorCount === 0) return null;

  return (
    <button
      aria-label={`${errorCount} form errors`}
      onClick={onClick}
      className="fixed bottom-4 right-4 p-3 rounded-full shadow-lg focus:outline-none"
    >
      <FaExclamationCircle size={24} aria-hidden="true" />
      <span className="sr-only">{errorCount} errors</span>
      <span className="absolute top-0 right-0 inline-flex items-center justify-center px-2 py-1 text-xs font-bold rounded-full bg-red-600 text-white">
        {errorCount}
      </span>
    </button>
  );
}
