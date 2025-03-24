import React, { useRef, useEffect } from "react";

interface ErrorDrawerProps {
  isOpen: boolean;
  errorsByTab: { index: number; count: number; message: string }[];
  onSelectTab: (i: number) => void;
  onClose: () => void;
}

export default function ErrorDrawer({
  isOpen,
  errorsByTab,
  onSelectTab,
  onClose,
}: ErrorDrawerProps) {
  // ⚠️ Specify the element type and initial value
  const drawerRef = useRef<HTMLUListElement | null>(null);

  useEffect(() => {
    if (isOpen && drawerRef.current) {
      drawerRef.current.focus();
    }
  }, [isOpen]);

  if (!isOpen) return null;

  return (
    <div
      role="dialog"
      aria-modal="true"
      className="fixed inset-0 bg-black bg-opacity-50 flex justify-center items-end"
      onClick={onClose}
    >
      <ul
        ref={drawerRef}
        className="bg-white w-full max-h-1/2 overflow-auto p-4 rounded-t-lg"
        tabIndex={-1}
        onClick={(e) => e.stopPropagation()}
      >
        {errorsByTab.map(({ index, count, message }) => (
          <li key={index}>
            <button
              onClick={() => { onSelectTab(index); onClose(); }}
              className="w-full text-left py-2"
            >
              Tab {index + 1}: {count} error{count > 1 && "s"} — {message}
            </button>
          </li>
        ))}
      </ul>
    </div>
  );
}
