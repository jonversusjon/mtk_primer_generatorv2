import React from "react";
import "../../styles/RestrictionSiteSummary.css";

function RestrictionSiteSummary({ sites }) {
  if (!sites || sites.length === 0) return null;

  return (
    <div className="restriction-sites-summary">
      <h3>Restriction Site Summary</h3>
      <table>
        <thead>
          <tr>
            <th>Enzyme</th>
            <th>Sequence</th>
            <th>Position</th>
            <th>Strand</th>
          </tr>
        </thead>
        <tbody>
          {sites.map((site, index) => (
            <tr key={index}>
              <td>{site.enzyme}</td>
              <td>{site.sequence}</td>
              <td>{site.position}</td>
              <td>{site.strand}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

export default RestrictionSiteSummary;
