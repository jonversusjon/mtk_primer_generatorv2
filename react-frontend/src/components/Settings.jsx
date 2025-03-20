import React, { useRef } from "react";
import "../styles/Settings.css";

const sliderValues = ["one", "a few", "many", "most", "all"];

// TODO: Add PCR reaction settings
/*
tm_threshold: 55.0
This likely refers to a threshold for melting temperature (Tm). Some settings
in Primer3 allow constraints based on a minimum threshold for Tm to ensure
primer stability.

min_3p_match: 10
This specifies the minimum number of consecutive base matches at the 3' end
between a primer and its target sequence. A higher value increases specificity,
reducing the chance of non-specific binding.

max_mismatches: 1 *** be sure this is tied to max_mut_per_site so it is never
                      less than max_mut_per_site
The maximum number of allowed mismatches between the primer and the target
sequence. A lower number increases primer specificity.

mv_conc: 50.0
The monovalent cation concentration (in mM), typically sodium (Na⁺) or
potassium (K⁺). This affects DNA stability and primer Tm calculations. The
default for PCR is often 50 mM Na⁺.

dv_conc: 1.5
The divalent cation concentration (in mM), typically magnesium (Mg²⁺). This
is crucial for DNA polymerase activity and affects primer binding. The default
Mg²⁺ concentration in PCR is usually 1.5–2.5 mM.

dntp_conc: 0.2
The dNTP (deoxynucleotide triphosphate) concentration (in mM). This value refers
to the concentration of each dNTP (dATP, dTTP, dCTP, dGTP) in the reaction. A
typical PCR concentration is 0.2 mM per dNTP.

dna_conc: 250.0
The DNA template concentration (in nM). This represents the assumed concentration
of single-stranded DNA in the PCR reaction, which impacts Tm calculations.

min_tm: 57
The minimum acceptable melting temperature (Tm) for a primer. This ensures that
primers are stable and will efficiently anneal under PCR conditions.
*/
function Settings({
  show,
  onClose,
  formData,
  updateFormData,
  availableSpecies,
}) {
  // Create a ref for the modal content - hooks must be called unconditionally
  const modalContentRef = useRef(null);

  // Handle clicks outside the modal content
  const handleOutsideClick = (e) => {
    // If we click outside the modal content, close the modal
    if (
      modalContentRef.current &&
      !modalContentRef.current.contains(e.target)
    ) {
      onClose();
    }
  };

  // If not shown, return null - but AFTER calling hooks
  if (!show) return null;

  return (
    <div className="settings-modal" onClick={handleOutsideClick}>
      <div className="settings-content" ref={modalContentRef}>
        {/* Species Selection */}
        <div className="form-group">
          <label htmlFor="species-select">Species:</label>
          <select
            id="species-select"
            className="form-control"
            value={formData.species}
            onChange={(e) => updateFormData("species", e.target.value)}
          >
            {availableSpecies.map((species) => (
              <option key={species} value={species}>
                {species}
              </option>
            ))}
          </select>
        </div>

        {/* Kozak Selection */}
        <div className="form-group">
          <label htmlFor="kozak-select">Kozak:</label>
          <select
            id="kozak-select"
            className="form-control"
            value={formData.kozak}
            onChange={(e) => updateFormData("kozak", e.target.value)}
          >
            <option value="MTK">MTK</option>
            <option value="Canonical">Canonical</option>
          </select>
        </div>

        {/* Mutations Setting */}
        <div className="form-group">
          <label htmlFor="mutations-select">Max mutations per site:</label>
          <select
            id="mutations-select"
            className="form-control"
            value={formData.max_mut_per_site}
            onChange={(e) =>
              updateFormData("max_mut_per_site", parseInt(e.target.value))
            }
          >
            <option value="1">1</option>
            <option value="2">2</option>
            <option value="3">3</option>
          </select>
        </div>

        {/* Results Limit Setting (Replaced Select with Slider) */}
        <div className="form-group">
          <label htmlFor="results-slider">Number of Results:</label>
          <input
            type="range"
            id="results-slider"
            className="form-control"
            min="0"
            max="4"
            step="1"
            value={sliderValues.indexOf(formData.results_limit)}
            onChange={(e) =>
              updateFormData(
                "results_limit",
                sliderValues[parseInt(e.target.value)]
              )
            }
          />
          <span style={{ marginLeft: "10px" }}>{formData.results_limit}</span>
        </div>

        {/* Verbose Mode Toggle */}
        <div className="form-check">
          <input
            type="checkbox"
            className="form-check-input"
            id="verbose-mode"
            checked={formData.verbose_mode}
            onChange={(e) => updateFormData("verbose_mode", e.target.checked)}
          />
          <label className="form-check-label" htmlFor="verbose-mode">
            Verbose
          </label>
        </div>

        {/* No footer with close button - will close by clicking outside or toggling settings button */}
      </div>
    </div>
  );
}

export default Settings;
