import React, { useState, useEffect } from "react";
import "../../styles/ProtocolTracker.css";
import RestrictionSiteSummary from "./RestrictionSiteSummary";

const ProgressStep = ({ name, progress, message }) => (
  <div className="border rounded-lg p-4 mb-4 shadow-sm bg-blue-50 border-blue-200">
    <div className="flex justify-between items-center mb-2">
      <span className="font-semibold text-gray-800">{name}</span>
      <span className="text-sm text-gray-500">{progress}%</span>
    </div>
    <div className="w-full bg-gray-200 rounded-full h-2 mb-2">
      <div
        className="bg-blue-500 h-2 rounded-full transition-all duration-300"
        style={{ width: `${progress}%` }}
      ></div>
    </div>
    {message && <p className="text-sm text-gray-600 mt-1">{message}</p>}
  </div>
);

const WaitingStep = ({ name }) => (
  <div className="border rounded-lg p-4 mb-4 shadow-sm bg-white border-gray-200">
    <span className="font-semibold text-gray-800">{name}</span>
  </div>
);

const TabButton = ({ name, isActive, onClick, notificationCount }) => (
  <button
    className={`protocol-tab-button ${isActive ? "active" : ""}`}
    onClick={onClick}
  >
    <div className="tab-icon-container">
      {notificationCount > 0 ? (
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
  <div className="bg-blue-50 border border-blue-200 rounded-md p-3 mb-4">
    <div className="flex items-start">
      <div className="flex-shrink-0 pt-0.5">
        <svg
          className="h-5 w-5 text-blue-500"
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
        <p className="text-sm text-blue-700">{message}</p>
      </div>
    </div>
  </div>
);

const TabContent = ({ stepName, stepData, messages, activeStep, sseData }) => {
  const stepMessages = messages.filter((msg) => msg.startsWith(`${stepName}:`));
  const [isMessagesOpen, setIsMessagesOpen] = useState(false);

  // Check if there's a display message for this step in the SSE data
  const displayMessage =
    sseData && sseData.step === stepName && sseData.display_message
      ? sseData.display_message
      : null;

  const renderStepSpecificContent = () => {
    // If there's a display message for this step, show it first
    if (displayMessage) {
      return (
        <div className="p-4">
          <DisplayMessage message={displayMessage} />
          {renderStepDetailContent()}
        </div>
      );
    }

    // Otherwise, just show regular step content
    return <div className="p-4">{renderStepDetailContent()}</div>;
  };

  const renderStepDetailContent = () => {
    if (
      stepName === "Restriction Site Detection" &&
      stepData?.restrictionSites?.length > 0
    ) {
      // Use the RestrictionSiteSummary component
      return <RestrictionSiteSummary sites={stepData.restrictionSites} />;
    }

  };

  return (
    <div className="protocol-tab-content">
      {activeStep && activeStep.name === stepName && (
        <ProgressStep
          name={activeStep.name}
          progress={activeStep.progress}
          message={activeStep.message}
        />
      )}
      {renderStepSpecificContent()}

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
                  <div key={index} className="text-gray-600 text-xs">
                    {msg.replace(`${stepName}: `, "")}
                  </div>
                ))}
              </div>
              <button
                className="text-blue-500 hover:text-blue-700 text-xs cursor-pointer mt-2"
                onClick={() => setIsMessagesOpen(false)}
                style={{
                  textDecoration: "none",
                  color: "inherit",
                  fontStyle: "italic",
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

  const [activeTab, setActiveTab] = useState(
    completedSteps.length > 0
      ? completedSteps[completedSteps.length - 1].name
      : activeSteps.length > 0
      ? activeSteps[0].name
      : null
  );

  // Set the active tab to the step that has an SSE message when it arrives
  // BUT only if one doesn't exist already - to prevent flashing
  useEffect(() => {
    if (sseData && sseData.step && sseData.display_message) {
      // Only update if we don't have an active tab yet or
      // if the SSE message is for the currently active tab
      if (!activeTab || sseData.step === activeTab) {
        const relevantStep = [...completedSteps, ...activeSteps].find(
          (step) => step.name === sseData.step
        );

        if (relevantStep) {
          setActiveTab(relevantStep.name);
        }
      }
    }
  }, [sseData, completedSteps, activeSteps, activeTab]);

  return (
    <div>
      <div className="protocol-tab-container">
        {completedSteps.concat(activeSteps).map((step) => (
          <TabButton
            key={step.name}
            name={step.name}
            isActive={activeTab === step.name}
            onClick={() => setActiveTab(step.name)}
            notificationCount={step.notificationCount}
          />
        ))}
      </div>

      {activeTab && (
        <TabContent
          stepName={activeTab}
          stepData={resultData[activeTab.replace(/\s+/g, "")]}
          messages={messages}
          activeStep={activeSteps.find((step) => step.name === activeTab)}
          sseData={sseData}
        />
      )}

      {waitingSteps.map((step) => (
        <WaitingStep key={step.name} name={step.name} />
      ))}
    </div>
  );
};

export default ProtocolTracker;
