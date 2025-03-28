// ProtocolTracker.jsx - Fixed to ensure steps are displayed as tabs
import React, { useState, useEffect } from "react";
import "../../styles/ProtocolTracker.css";
import RestrictionSiteSummary from "./RestrictionSiteSummary";
import MutationAnalysisSummary from "./MutationAnalysisSummary";

// UI Helper Components
const ProgressStep = ({ name, stepProgress, message }) => (
  <div className="border rounded-lg p-4 mb-4 shadow-sm bg-blue-50 dark:bg-blue-900/20 border-blue-200 dark:border-blue-800">
    <div className="flex justify-between items-center mb-2">
      <span className="font-semibold text-gray-800 dark:text-gray-100">
        {name}
      </span>
      <span className="text-sm text-gray-500 dark:text-gray-400">
        {stepProgress}%
      </span>
    </div>
    <div className="w-full bg-gray-200 dark:bg-gray-700 rounded-full h-2 mb-2">
      <div
        className="bg-blue-500 dark:bg-blue-400 h-2 rounded-full transition-all duration-300"
        style={{ width: `${stepProgress}%` }}
      ></div>
    </div>
    {message && (
      <p className="text-sm text-gray-600 dark:text-gray-300 mt-1">{message}</p>
    )}
  </div>
);

// const WaitingStep = ({ name }) => (
//   <div className="border rounded-lg p-4 mb-4 shadow-sm bg-white dark:bg-gray-800 border-gray-200 dark:border-gray-700">
//     <span className="font-semibold text-gray-800 dark:text-gray-200">
//       {name}
//     </span>
//   </div>
// );

// Tab Button Component - This renders the actual tab buttons
const TabButton = ({
  name,
  isActive,
  onClick,
  notificationCount,
  status,
  stepProgress,
}) => (
  <button
    className={`protocol-tab-button ${isActive ? "active" : ""}`}
    onClick={onClick}
  >
    {notificationCount > 0 ? (
      <span className="tab-notification">{notificationCount}</span>
    ) : status === "completed" || stepProgress === 100 ? (
      <div className="checkmark-container">✓</div>
    ) : (
      <div className="checkmark-placeholder"></div>
    )}
    <span className="tab-label">{name}</span>
  </button>
);

const DisplayMessage = ({ message, timestamp }) => (
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
        <div>
          <p className="text-sm text-blue-700 dark:text-blue-300">{message}</p>
          {timestamp && (
            <p className="text-xs text-gray-500 dark:text-gray-400 mt-1">
              {new Date(timestamp).toLocaleTimeString()}
            </p>
          )}
        </div>
      </div>
    </div>
  </div>
);

// Tab Content Component
const TabContent = ({ stepName, messages, activeStep, sseData, callouts }) => {
  // Filter messages relevant to this step
  const stepMessages = Array.isArray(messages)
    ? messages.filter((m) => m.step === stepName).map((m) => m.message)
    : [];

  const [isMessagesOpen, setIsMessagesOpen] = useState(false);
  const [isPayloadVisible, setIsPayloadVisible] = useState(false);
  // es-lint-disable-next-line no-unused-vars
  // const [selectedMutationSetIndex, _setSelectedMutationSetIndex] = useState(0);

  const [selectedMutationSetIndex] = useState(0);

  // Get the SSE data for the current step
  const stepSseData = sseData?.[stepName];

  const renderStepDetailContent = () => {
    if (!stepSseData) {
      console.debug(`No SSE data found for step: ${stepName}`);
      return (
        <div className="p-4 text-gray-600 dark:text-gray-400 text-center">
          Waiting for data...
        </div>
      );
    }

    // Render content based on the step
    switch (stepName) {
      case "Restriction Sites":
        if (stepSseData.sites && stepSseData.sites.length > 0) {
          return <RestrictionSiteSummary sites={stepSseData.sites} />;
        }
        break;

      case "Mutation Analysis":
        // First check for detailed mutation analysis
        if (
          stepSseData.mutationSets &&
          stepSseData.mutationSets.length > 0 &&
          stepSseData.restrictionSites &&
          stepSseData.restrictionSites.length > 0
        ) {
          const restrictionSite = stepSseData.restrictionSites.find(
            (site) => site.context_seq && site.context_rs_indices && site.codons
          );
          if (restrictionSite) {
            return (
              <div className="mt-2">
                <h3 className="font-semibold text-gray-700 dark:text-gray-200 mb-4">
                  Mutation Analysis Summary
                </h3>
                <MutationAnalysisSummary
                  restrictionSite={restrictionSite}
                  mutationSets={stepSseData.mutationSets}
                  selectedMutationSetIndex={selectedMutationSetIndex}
                />
              </div>
            );
          }
        }

        // Fallback to basic mutations table
        if (stepSseData.mutations && stepSseData.mutations.length > 0) {
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
                    {stepSseData.mutations.map((mutation, idx) => (
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
        }
        break;

      case "Primer Design":
        if (
          (stepSseData.edgePrimers &&
            Object.keys(stepSseData.edgePrimers).length > 0) ||
          (stepSseData.mutPrimers &&
            Object.keys(stepSseData.mutPrimers).length > 0)
        ) {
          return (
            <div className="mt-2">
              <h3 className="font-semibold text-gray-700 dark:text-gray-200 mb-2">
                Designed Primers:
              </h3>
              {stepSseData.edgePrimers &&
                Object.keys(stepSseData.edgePrimers).length > 0 && (
                  <div className="mb-4">
                    <h4 className="text-sm font-medium text-gray-600 dark:text-gray-300 mb-1">
                      Edge Primers:
                    </h4>
                    <div className="border dark:border-gray-700 rounded overflow-hidden">
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
                          {Object.entries(stepSseData.edgePrimers).map(
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
              {stepSseData.mutPrimers &&
                Object.keys(stepSseData.mutPrimers).length > 0 && (
                  <div>
                    <h4 className="text-sm font-medium text-gray-600 dark:text-gray-300 mb-1">
                      Mutation Primers:
                    </h4>
                    <div className="border dark:border-gray-700 rounded overflow-hidden">
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
                          {Object.entries(stepSseData.mutPrimers).map(
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
        }
        break;

      case "PCR Reaction Grouping":
        if (stepSseData.pcrReactions && stepSseData.pcrReactions.length > 0) {
          return (
            <div className="mt-2">
              <h3 className="font-semibold text-gray-700 dark:text-gray-200 mb-2">
                PCR Reactions:
              </h3>
              <div className="space-y-4">
                {stepSseData.pcrReactions.map((reaction, idx) => (
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
        break;

      default:
        return (
          <div className="p-4 text-gray-600 dark:text-gray-400 text-center">
            Waiting for data...
          </div>
        );
    }

    return (
      <div className="p-4 text-gray-600 dark:text-gray-400 text-center">
        Processing data...
      </div>
    );
  };

  return (
    <div className="protocol-tab-content">
      {/* Show progress for the active step */}
      {activeStep && activeStep.name === stepName && (
        <ProgressStep
          name={activeStep.name}
          stepProgress={activeStep.stepProgress}
          message={activeStep.message}
        />
      )}

      <div className="p-4">
        {/* Display callouts if any */}
        {callouts && callouts.length > 0 && (
          <div className="space-y-2">
            {callouts
              .slice()
              .reverse()
              .map((callout, idx) => (
                <DisplayMessage
                  key={idx}
                  message={callout.message}
                  timestamp={callout.timestamp}
                />
              ))}
          </div>
        )}

        {/* Main content for this step */}
        {renderStepDetailContent()}

        {/* Raw SSE payload toggle */}
        {stepSseData && Object.keys(stepSseData).length > 0 && (
          <div className="border-t border-gray-200 mt-4 pt-4">
            <button
              className="text-blue-500 hover:text-blue-700 text-xs cursor-pointer subtle-link mb-2 dark:text-gray-400 dark:hover:text-gray-300"
              onClick={() => setIsPayloadVisible(!isPayloadVisible)}
              style={{
                textDecoration: "none",
                fontStyle: "italic",
                border: "none",
                background: "none",
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

      {/* Toggle display of all step-specific messages */}
      {stepMessages.length > 0 && (
        <div className="border-t border-gray-200 p-4">
          {!isMessagesOpen && (
            <button
              className="text-blue-500 hover:text-blue-700 text-xs cursor-pointer subtle-link"
              onClick={() => setIsMessagesOpen(true)}
              style={{
                textDecoration: "none",
                fontStyle: "italic",
                border: "none",
                background: "none",
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
                  fontStyle: "italic",
                  border: "none",
                  background: "none",
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

// Main ProtocolTracker Component
const ProtocolTracker = ({ steps, messages, sseData }) => {
  // Filter steps by status
  const completedSteps = steps.filter((step) => step.status === "completed");
  const activeSteps = steps.filter((step) => step.status === "active");
  // es-lint-disable-next-line no-unused-vars
  // const waitingSteps = steps.filter((step) => step.status === "waiting");

  // Flag for manual tab selection
  const [userSelectedTab, setUserSelectedTab] = useState(false);

  // Accumulated callouts per step
  const [stepCallouts, setStepCallouts] = useState({});

  // Determine the initial active tab
  const [activeTab, setActiveTab] = useState(() => {
    if (activeSteps.length > 0) return activeSteps[0].name;
    if (completedSteps.length > 0)
      return completedSteps[completedSteps.length - 1].name;
    return steps.length > 0 ? steps[0].name : null;
  });

  // Accumulate callouts from SSE data
  useEffect(() => {
    if (!sseData) return;

    Object.entries(sseData)
      .filter(([_, data]) => data?.callout)
      .forEach(([stepName, stepData]) => {
        setStepCallouts((prev) => {
          const existingCallouts = prev[stepName] || [];
          if (existingCallouts.some((c) => c.message === stepData.callout))
            return prev;
          return {
            ...prev,
            [stepName]: [
              ...existingCallouts,
              {
                message: stepData.callout,
                timestamp: stepData.timestamp || Date.now(),
              },
            ],
          };
        });
      });
  }, [sseData]);

  const handleTabSelect = (tabName) => {
    setActiveTab(tabName);
    setUserSelectedTab(true);
  };

  // Auto-switch tab if a new active step is detected and the user hasn't manually selected a tab
  useEffect(() => {
    if (activeSteps.length > 0 && !userSelectedTab) {
      const mostRecentActiveStep = activeSteps[activeSteps.length - 1].name;
      setActiveTab(mostRecentActiveStep);
    }
  }, [activeSteps, userSelectedTab]);

  // Auto-switch based on callouts if the user hasn't manually selected a tab
  useEffect(() => {
    if (userSelectedTab || !sseData) return;

    const stepsWithCallouts = Object.entries(sseData)
      .filter(([_, data]) => data?.callout)
      .map(([stepName, data]) => ({
        stepName,
        timestamp: data.timestamp || Date.now(),
      }))
      .sort((a, b) => b.timestamp - a.timestamp);

    if (stepsWithCallouts.length > 0) {
      const mostRecentStep = stepsWithCallouts[0].stepName;
      const isStepRendered = [...completedSteps, ...activeSteps].some(
        (step) => step.name === mostRecentStep
      );

      if (isStepRendered) {
        setActiveTab(mostRecentStep);
      }
    }
  }, [sseData, completedSteps, activeSteps, userSelectedTab]);

  // Only show tabs that are either completed or active
  const tabsToShow = steps.filter(
    (step) => step.status === "completed" || step.status === "active"
  );

  // Show at least 1 step, even if it's waiting
  const displaySteps = tabsToShow.length > 0 ? tabsToShow : steps.slice(0, 1);

  // If no steps are available at all, show a loading indicator
  if (steps.length === 0) {
    return (
      <div className="p-4 text-center text-gray-500">
        Loading protocol steps...
      </div>
    );
  }

  return (
    <div>
      {/* Tab Buttons Container - This is the tab navigation bar */}
      <div className="protocol-tab-container">
        {displaySteps.map((step) => (
          <TabButton
            key={step.name}
            name={step.name}
            isActive={activeTab === step.name}
            onClick={() => handleTabSelect(step.name)}
            notificationCount={step.notificationCount || 0}
            status={step.status}
            stepProgress={step.stepProgress || 0}
          />
        ))}
      </div>

      {/* Only render the active tab content */}
      {activeTab && (
        <div className="tab-content">
          <TabContent
            stepName={activeTab}
            messages={messages || []}
            activeStep={activeSteps.find((step) => step.name === activeTab)}
            sseData={sseData || {}}
            callouts={stepCallouts[activeTab] || []}
          />
        </div>
      )}
    </div>
  );
};

export default ProtocolTracker;
