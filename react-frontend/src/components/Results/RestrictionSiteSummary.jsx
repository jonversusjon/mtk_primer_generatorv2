import React from 'react';
import '../../styles/RestrictionSiteSummary.css';

function RestrictionSiteSummary({ sites }) {
  if (!sites || sites.length === 0) return null;
  
  return (
    <div className="restriction-site-summary">
      <h3 className="restriction-site-summary-table-title">Restriction Site Summary</h3>
      <div className="restriction-site-summary-table-container">
        <table>
          <thead>
            <tr>
              <th>Site Type</th>
              <th>Number of Instances</th>
              <th>Position(s)</th>
            </tr>
          </thead>
          <tbody>
            {sites.map((site, index) => (
              <tr key={index}>
                <td>{site.enzyme} Restriction Site</td>
                <td>{site.count}</td>
                <td>{site.positions.join(', ')}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

export default RestrictionSiteSummary;