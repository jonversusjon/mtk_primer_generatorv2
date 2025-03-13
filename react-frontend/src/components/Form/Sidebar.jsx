import React from "react";
import { Tooltip } from "react-tooltip";
import "../../styles/global/theme.css";
import "../../styles/Sidebar.css";

const Sidebar = ({
  sequences,
  errorsBySequence,
  onSelectTab,
  activeTabIndex,
}) => {
  return (
    <div className="sidebar">
      <h3>Sequences Overview</h3>
      <ul>
        {sequences.map((seq, index) => {
          const errors = errorsBySequence[index] || [];
          return (
            <li
              key={index}
              className={`sidebar-item ${
                activeTabIndex === index ? "active" : ""
              }`}
              onClick={() => onSelectTab(index)}
            >
              <div className="sidebar-item-info">
                <div className="primer-name">
                  {seq.primerName || "Unnamed Primer"}
                </div>
                <div className="sequence-length">
                  {(seq.sequence || "").length} bp
                </div>
              </div>

              {errors.length > 0 ? (
                <>
                  <div
                    data-tooltip-id={`error-tooltip-${index}`}
                    data-tooltip-content={errors.join("\n")}
                    className="error-badge-container"
                  >
                    <span className="error-badge">{errors.length}</span>
                  </div>
                  <Tooltip
                    id={`error-tooltip-${index}`}
                    place="right"
                    style={{ whiteSpace: "pre-line" }}
                  />
                </>
              ) : (
                <div className="checkmark-container">✔</div>
              )}
            </li>
          );
        })}
      </ul>
    </div>
  );
};

export default Sidebar;
