import React, { useState, useEffect } from "react";
import "../../styles/ProtocolTracker.css";
import RestrictionSiteSummary from "./RestrictionSiteSummary";

const ProgressStep = ({ name, progress, message }) => (
  <div className="border rounded-lg p-4 mb-4 shadow-sm bg-blue-50 dark:bg-blue-900/20 border-blue-200 dark:border-blue-800">
    <div className="flex justify-between items-center mb-2">
      <span className="font-semibold text-gray-800 dark:text-gray-100">
        {name}
      </span>
      <span className="text-sm text-gray-500 dark:text-gray-400">
        {progress}%
      </span>
    </div>
    <div className="w-full bg-gray-200 dark:bg-gray-700 rounded-full h-2 mb-2">
      <div
        className="bg-blue-500 dark:bg-blue-400 h-2 rounded-full transition-all duration-300"
        style={{ width: `${progress}%` }}
      ></div>
    </div>
    {message && (
      <p className="text-sm text-gray-600 dark:text-gray-300 mt-1">{message}</p>
    )}
  </div>
);

const WaitingStep = ({ name }) => (
  <div className="border rounded-lg p-4 mb-4 shadow-sm bg-white dark:bg-gray-800 border-gray-200 dark:border-gray-700">
    <span className="font-semibold text-gray-800 dark:text-gray-200">
      {name}
    </span>
  </div>
);

const TabButton = ({ name, isActive, onClick, notificationCount }) => (
  <button
    className={`protocol-tab-button ${isActive ? "active" : ""}`}
    onClick={onClick}
  >
    <div className="tab-icon-container">
      {/* Display checkmark only if completed (progress 100) and no notifications */}
      {notificationCount === 0 &&
      isActive /* Assuming isActive implies it could be completed or active */ ? (
        <div className="checkmark-container">
          ✔
        </div> /* Or refine based on actual step status if available */
      ) : notificationCount > 0 ? (
        <span className="tab-notification">{notificationCount}</span>
      ) : (
        <div className="checkmark-container">✔</div>
      )}
    </div>
    <span className="tab-label">{name}</span>
  </button>
);

// DisplayMessage component for showing SSE display_messages
const DisplayMessage = ({ message }) => (
  <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded-md p-3 mb-4">
    <div className="flex items-start">
      <div className="flex-shrink-0 pt-0.5">
        <svg
          className="h-5 w-5 text-blue-500 dark:text-blue-400"
          xmlns="http://www.w3.org/2000/svg"
          viewBox="0 0 20 20"
          fill="currentColor"
        >
          <path
            fillRule="evenodd"
            d="M18 10a8 8 0 11-16 0 8 8 0 0116 0zm-7-4a1 1 0 11-2 0 1 1 0 012 0zM9 9a1 1 0 000 2v3a1 1 0 001 1h1a1 1 0 100-2h-1V9z"
            clipRule="evenodd"
          />
        </svg>
      </div>
      <div className="ml-3">
        <p className="text-sm text-blue-700 dark:text-blue-300">{message}</p>
      </div>
    </div>
  </div>
);

// Map of step names to resultData keys
const stepToDataKeyMap = {
  Preprocessing: "Preprocessing",
  "Restriction Site Detection": "RestrictionSiteDetection",
  "Mutation Analysis": "MutationAnalysis",
  "Primer Design": "PrimerDesign",
  "PCR Reaction Grouping": "PCRReactionGrouping",
};

const TabContent = ({ stepName, stepData, messages, activeStep, sseData }) => {
  const stepMessages = messages.filter((msg) => msg.startsWith(`${stepName}:`));
  const [isMessagesOpen, setIsMessagesOpen] = useState(false);
  const [isPayloadVisible, setIsPayloadVisible] = useState(false); // State for payload visibility

  // Get step-specific SSE data instead of global sseData
  const stepSseData = sseData ? sseData[stepName] : null;

  // Determine if there's a specific callout message for this step
  const callout =
    stepSseData && stepSseData.callout ? stepSseData.callout : null;

  // Helper function to get the specific content *element* for the current step
  const renderStepDetailContent = () => {
    // Get the appropriate data key for this step
    const dataKey = stepToDataKeyMap[stepName];
    if (!dataKey || !stepData || !stepData[dataKey]) {
      console.warn(`No data found for step: ${stepName}, key: ${dataKey}`);
      return null;
    }

    // Get data specific to this step
    const data = stepData[dataKey];

    if (stepName === "Restriction Sites" && data.restrictionSites?.length > 0) {
      return <RestrictionSiteSummary sites={data.restrictionSites} />;
    }
    // Add other 'else if' conditions here for different step names
    else if (stepName === "Mutation Analysis" && data.mutations?.length > 0) {
      // You would need to create a component for displaying mutations
      // return <MutationAnalysisSummary mutations={data.mutations} />;
      return (
        <div className="mt-2">
          <h3 className="font-semibold text-gray-700 dark:text-gray-200 mb-2">
            Mutations Found:
          </h3>
          <div className="border dark:border-gray-700 rounded overflow-hidden">
            <table className="min-w-full divide-y divide-gray-200">
              <thead className="bg-gray-50 dark:bg-gray-700">
                <tr>
                  <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 dark:text-gray-300 uppercase tracking-wider">
                    Type
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 dark:text-gray-300 uppercase tracking-wider">
                    Position
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 dark:text-gray-300 uppercase tracking-wider">
                    Sequence
                  </th>
                </tr>
              </thead>
              <tbody className="bg-white dark:bg-gray-800 divide-y divide-gray-200 dark:divide-gray-700">
                {data.mutations.map((mutation, idx) => (
                  <tr key={idx}>
                    <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500 dark:text-gray-300">
                      {mutation.type}
                    </td>
                    <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500 dark:text-gray-300">
                      {mutation.position}
                    </td>
                    <td className="px-6 py-4 whitespace-nowrap text-sm font-mono text-gray-500 dark:text-gray-300">
                      {mutation.sequence}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      );
    } else if (
      stepName === "Primer Design" &&
      (data.edgePrimers || data.mutPrimers)
    ) {
      // Placeholder for primer design summary
      return (
        <div className="mt-2">
          <h3 className="font-semibold text-gray-700 mb-2">
            Designed Primers:
          </h3>
          {data.edgePrimers && Object.keys(data.edgePrimers).length > 0 && (
            <div className="mb-4">
              <h4 className="text-sm font-medium text-gray-600 mb-1">
                Edge Primers:
              </h4>
              <div className="border rounded overflow-hidden">
                <table className="min-w-full divide-y divide-gray-200">
                  <thead className="bg-gray-50 dark:bg-gray-700">
                    <tr>
                      <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 dark:text-gray-300 uppercase tracking-wider">
                        Name
                      </th>
                      <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 dark:text-gray-300 uppercase tracking-wider">
                        Sequence
                      </th>
                    </tr>
                  </thead>
                  <tbody className="bg-white dark:bg-gray-800 divide-y divide-gray-200 dark:divide-gray-700">
                    {Object.entries(data.edgePrimers).map(
                      ([name, primer], idx) => (
                        <tr key={idx}>
                          <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500 dark:text-gray-300">
                            {name}
                          </td>
                          <td className="px-6 py-4 whitespace-nowrap text-sm font-mono text-gray-500 dark:text-gray-300">
                            {primer.sequence || primer}
                          </td>
                        </tr>
                      )
                    )}
                  </tbody>
                </table>
              </div>
            </div>
          )}
          {data.mutPrimers && Object.keys(data.mutPrimers).length > 0 && (
            <div>
              <h4 className="text-sm font-medium text-gray-600 mb-1">
                Mutation Primers:
              </h4>
              <div className="border rounded overflow-hidden">
                <table className="min-w-full divide-y divide-gray-200">
                  <thead className="bg-gray-50 dark:bg-gray-700">
                    <tr>
                      <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 dark:text-gray-300 uppercase tracking-wider">
                        Name
                      </th>
                      <th className="px-6 py-3 text-left text-xs font-medium text-gray-500 dark:text-gray-300 uppercase tracking-wider">
                        Sequence
                      </th>
                    </tr>
                  </thead>
                  <tbody className="bg-white dark:bg-gray-800 divide-y divide-gray-200 dark:divide-gray-700">
                    {Object.entries(data.mutPrimers).map(
                      ([name, primer], idx) => (
                        <tr key={idx}>
                          <td className="px-6 py-4 whitespace-nowrap text-sm text-gray-500 dark:text-gray-300">
                            {name}
                          </td>
                          <td className="px-6 py-4 whitespace-nowrap text-sm font-mono text-gray-500 dark:text-gray-300">
                            {primer.sequence || primer}
                          </td>
                        </tr>
                      )
                    )}
                  </tbody>
                </table>
              </div>
            </div>
          )}
        </div>
      );
    } else if (
      stepName === "PCR Reaction Grouping" &&
      data.pcrReactions?.length > 0
    ) {
      // Placeholder for PCR reactions display
      return (
        <div className="mt-2">
          <h3 className="font-semibold text-gray-700 mb-2">PCR Reactions:</h3>
          <div className="space-y-4">
            {data.pcrReactions.map((reaction, idx) => (
              <div
                key={idx}
                className="border rounded p-3 bg-gray-50 dark:bg-gray-800 dark:border-gray-700"
              >
                <h4 className="font-medium text-gray-700 dark:text-gray-200 mb-2">
                  Reaction {idx + 1}
                </h4>
                <div className="grid grid-cols-2 gap-2 text-sm">
                  <div className="col-span-2">
                    <span className="font-medium dark:text-gray-300">
                      Template:
                    </span>{" "}
                    <span className="dark:text-gray-300">
                      {reaction.template || "N/A"}
                    </span>
                  </div>
                  <div>
                    <span className="font-medium dark:text-gray-300">
                      Forward Primer:
                    </span>{" "}
                    <span className="dark:text-gray-300">
                      {reaction.forwardPrimer || "N/A"}
                    </span>
                  </div>
                  <div>
                    <span className="font-medium dark:text-gray-300">
                      Reverse Primer:
                    </span>{" "}
                    <span className="dark:text-gray-300">
                      {reaction.reversePrimer || "N/A"}
                    </span>
                  </div>
                  {reaction.product && (
                    <div className="col-span-2">
                      <span className="font-medium dark:text-gray-300">
                        Product:
                      </span>
                      <div className="font-mono text-xs mt-1 p-1 bg-gray-100 dark:bg-gray-900 rounded dark:text-gray-300">
                        {reaction.product}
                      </div>
                    </div>
                  )}
                </div>
              </div>
            ))}
          </div>
        </div>
      );
    }

    // Return null or a placeholder if no specific content for this step
    return null;
  };

  // Get the specific detail content using the helper
  const detailContent = renderStepDetailContent();

  return (
    <div className="protocol-tab-content">
      {/* Display Progress if this is the active step */}
      {activeStep && activeStep.name === stepName && (
        <ProgressStep
          name={activeStep.name}
          progress={activeStep.progress}
          message={activeStep.message}
        />
      )}

      {/* Main content area for the tab */}
      <div className="p-4">
        {/* Display Callout message if it exists */}
        {callout && <DisplayMessage message={callout} />}

        {/* Display the specific detail content */}
        {detailContent}

        {/* Section to display SSE Payload (Toggleable) */}
        {/* Only show payload section if stepSseData exists for this specific step */}
        {stepSseData && (
          <div className="border-t border-gray-200 mt-4 pt-4">
            <button
              className="text-blue-500 hover:text-blue-700 text-xs cursor-pointer subtle-link mb-2 dark:text-gray-400 dark:hover:text-gray-300"
              onClick={() => setIsPayloadVisible(!isPayloadVisible)}
              style={{
                textDecoration: "none",
                color: "inherit",
                fontStyle: "italic",
                border: "none",
              }}
            >
              {isPayloadVisible
                ? "Hide Raw SSE Payload"
                : "Show Raw SSE Payload"}
            </button>
            {isPayloadVisible && (
              <pre className="bg-gray-100 dark:bg-gray-800 p-2 rounded text-xs overflow-x-auto dark:text-gray-300">
                {JSON.stringify(stepSseData, null, 2)}
              </pre>
            )}
          </div>
        )}
      </div>

      {/* Section for displaying all step-specific messages (Toggleable) */}
      {stepMessages.length > 0 && (
        <div className="border-t border-gray-200 p-4">
          {!isMessagesOpen && (
            <button
              className="text-blue-500 hover:text-blue-700 text-xs cursor-pointer subtle-link"
              onClick={() => setIsMessagesOpen(true)}
              style={{
                textDecoration: "none",
                color: "inherit",
                fontStyle: "italic",
                border: "none",
              }}
            >
              see all messages
            </button>
          )}
          {isMessagesOpen && (
            <>
              <div className="space-y-1 max-h-40 overflow-y-auto">
                {stepMessages.map((msg, index) => (
                  <div
                    key={index}
                    className="text-gray-600 dark:text-gray-400 text-xs"
                  >
                    {msg.replace(`${stepName}: `, "")}
                  </div>
                ))}
              </div>
              <button
                className="text-blue-500 hover:text-blue-700 text-xs cursor-pointer mt-2 dark:text-gray-400 dark:hover:text-gray-300"
                onClick={() => setIsMessagesOpen(false)}
                style={{
                  textDecoration: "none",
                  color: "inherit",
                  fontStyle: "italic",
                  border: "none",
                }}
              >
                hide messages
              </button>
            </>
          )}
        </div>
      )}
    </div>
  );
};

const ProtocolTracker = ({ steps, messages, resultData, sseData }) => {
  const completedSteps = steps.filter((step) => step.status === "completed");
  const activeSteps = steps.filter((step) => step.status === "active");
  const waitingSteps = steps.filter((step) => step.status === "waiting");

  // Determine initial active tab: last completed or first active
  const [activeTab, setActiveTab] = useState(() => {
    if (activeSteps.length > 0) return activeSteps[0].name;
    if (completedSteps.length > 0)
      return completedSteps[completedSteps.length - 1].name;
    return steps.length > 0 ? steps[0].name : null; // Fallback to first step if none are active/completed
  });

  // Effect to potentially switch tab based on SSE *callout* message
  // (Consider if this behavior is desired - might jump user unexpectedly)
  useEffect(() => {
    // Only attempt to switch if sseData has new information
    if (sseData) {
      // Find the step that has the most recent SSE data with a callout
      const stepsWithCallouts = Object.entries(sseData)
        .filter(([_, data]) => data?.callout)
        .map(([stepName, data]) => ({
          stepName,
          timestamp: data.timestamp || Date.now(), // Use timestamp if available
        }))
        .sort((a, b) => b.timestamp - a.timestamp); // Sort by most recent

      if (stepsWithCallouts.length > 0) {
        // Get the most recent step with a callout
        const mostRecentStep = stepsWithCallouts[0].stepName;

        // Check if this step is currently rendered
        const isStepRendered = [...completedSteps, ...activeSteps].some(
          (step) => step.name === mostRecentStep
        );

        if (isStepRendered) {
          setActiveTab(mostRecentStep);
        }
      }
    }
  }, [sseData, completedSteps, activeSteps]);

  // Derive the list of tabs to show (completed + active)
  const tabsToShow = steps.filter(
    (step) => step.status === "completed" || step.status === "active"
  );

  return (
    <div>
      <div className="protocol-tab-container">
        {tabsToShow.map((step) => (
          <TabButton
            key={step.name}
            name={step.name}
            isActive={activeTab === step.name}
            onClick={() => setActiveTab(step.name)}
            notificationCount={step.notificationCount || 0}
          />
        ))}
      </div>

      {activeTab && (
        <TabContent
          stepName={activeTab}
          stepData={resultData} // Pass the full resultData object
          messages={messages}
          activeStep={activeSteps.find((step) => step.name === activeTab)}
          sseData={sseData} // Pass the entire sseData object
        />
      )}

      {waitingSteps.map((step) => (
        <WaitingStep key={step.name} name={step.name} />
      ))}
    </div>
  );
};

export default ProtocolTracker;
