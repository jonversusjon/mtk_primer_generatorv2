import React from "react";

function RestrictionSiteSummary({ restrictionSites }) {
  if (!restrictionSites || restrictionSites.length === 0) return null;

  return (
    <div className="restriction-sites-summary section-container">
      <h3>Internal BsaI/BsmbI sites found</h3>
      <div className="restriction-site-summary-table-wrapper">
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
            {restrictionSites.map((site, index) => (
              <tr key={index}>
                <td>{site.enzyme}</td>
                <td>{site.recognitionSeq}</td>
                <td>{site.position}</td>
                <td>{site.strand}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

export default RestrictionSiteSummary;
