import React from 'react';
import '../../styles/Settings.css';

function Settings({ show, onClose, formData, updateFormData, availableSpecies }) {
  if (!show) return null;
  
  return (
    <div id="settings-card" className="settings-card">
      <div className="settings-content">
        {/* Species Selection */}
        <div className="form-group">
          <label htmlFor="species">Species:</label>
          <select
            id="species"
            name="species"
            className="form-select"
            value={formData.species}
            onChange={(e) => updateFormData('species', e.target.value)}
          >
            <option value="">Select species...</option>
            {availableSpecies.map((species) => (
              <option key={species} value={species}>
                {species}
              </option>
            ))}
          </select>
        </div>

        {/* Kozak Selection */}
        <div className="form-group">
          <label htmlFor="kozak">Kozak:</label>
          <select
            id="kozak"
            name="kozak"
            className="form-select"
            value={formData.kozak}
            onChange={(e) => updateFormData('kozak', e.target.value)}
          >
            <option value="MTK">MTK</option>
            <option value="Canonical">Canonical</option>
          </select>
        </div>

        {/* Mutations Setting */}
        <div className="form-group">
          <label htmlFor="maxMutPerSite">Max mutations per site:</label>
          <select
            id="maxMutPerSite"
            name="max_mut_per_site"
            className="form-select"
            value={formData.max_mut_per_site}
            onChange={(e) => updateFormData('max_mut_per_site', parseInt(e.target.value))}
          >
            <option value="1">1</option>
            <option value="2">2</option>
            <option value="3">3</option>
          </select>
        </div>

        {/* Verbose Mode Toggle */}
        <div className="form-group checkbox">
          <input
            type="checkbox"
            id="verbose_mode"
            name="verbose_mode"
            checked={formData.verbose_mode}
            onChange={(e) => updateFormData('verbose_mode', e.target.checked)}
          />
          <label htmlFor="verbose_mode">Verbose</label>
        </div>
        
        <button 
          id="close-settings" 
          className="close-settings-btn" 
          onClick={onClose}
        >
          Close Settings
        </button>
      </div>
    </div>
  );
}

export default Settings;