// Results/RawSseDisplay.jsx
import React from 'react';

// Simple component to render a list of SSE events using <details> for toggling
const RawSseDisplay = ({ sseEvents }) => {

  // If there are no events, display a message
  if (!sseEvents || sseEvents.length === 0) {
    return <p className="text-gray-500 dark:text-gray-400">Waiting for Server-Sent Events...</p>;
  }

  return (
    <div className="space-y-2"> {/* Add vertical space between event details */}
      <h3 className="text-lg font-semibold mb-3 dark:text-gray-200">Raw Server-Sent Events</h3>
      {/* Map over the received SSE events */}
      {sseEvents.map((event, index) => (
        <details key={index} className="border rounded-md border-gray-300 dark:border-gray-600 overflow-hidden">
          {/* The summary part - visible when collapsed */}
          <summary className="cursor-pointer p-3 bg-gray-100 dark:bg-gray-700 hover:bg-gray-200 dark:hover:bg-gray-600 flex justify-between items-center">
            <span className="font-medium text-sm dark:text-gray-200">
              Event {index + 1}: {event.step || 'No Step Info'} {/* Show step name if available */}
            </span>
            <span className="text-xs text-gray-500 dark:text-gray-400">
              {/* Display timestamp if available */}
              {event.clientTimestamp ? new Date(event.clientTimestamp).toLocaleTimeString() : 'No Timestamp'}
            </span>
          </summary>
          {/* The content part - visible when expanded */}
          <div className="p-3 bg-white dark:bg-gray-800">
            <pre className="text-xs overflow-x-auto whitespace-pre-wrap break-all bg-gray-50 dark:bg-gray-900 p-2 rounded dark:text-gray-300">
              {/* Display the full event data nicely formatted as JSON */}
              {JSON.stringify(event, null, 2)}
            </pre>
          </div>
        </details>
      ))}
    </div>
  );
};

export default RawSseDisplay;