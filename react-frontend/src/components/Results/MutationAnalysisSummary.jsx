import React, { useState, useEffect } from 'react';

/**
 * MutationAnalysisSummary - A component for visualizing restriction site mutations
 * 
 * This component displays:
 * 1. Context sequence with highlighted recognition site
 * 2. Original codons within the recognition site
 * 3. Alternative codons with their usage frequencies
 */
const MutationAnalysisSummary = ({ 
  restrictionSite, 
  mutationSets, 
  selectedMutationSetIndex = 0 
}) => {
  const [selectedMutationSet, setSelectedMutationSet] = useState(null);

  useEffect(() => {
    if (mutationSets && mutationSets.length > 0) {
      setSelectedMutationSet(mutationSets[selectedMutationSetIndex]);
    }
  }, [mutationSets, selectedMutationSetIndex]);

  if (!restrictionSite || !selectedMutationSet) {
    return <div className="p-4">Loading mutation analysis data...</div>;
  }

  // Destructure the restriction site data
  const { 
    context_seq, 
    context_rs_indices, 
    context_first_base, 
    codons: originalCodons 
  } = restrictionSite;

  // Get mutations from the selected mutation set
  const mutations = selectedMutationSet.mutations || [];

  // Function to get alternative codons for a specific original codon position
  const getAlternativeCodon = (codonPosition) => {
    for (const mutation of mutations) {
      for (const mutCodon of (mutation.mut_codons || [])) {
        if (mutCodon.nth_codon_in_rs === codonPosition) {
          return mutCodon.codon;
        }
      }
    }
    return null;
  };

  // Prepare the context sequence with its styling
  const renderContextSequence = () => {
    return context_seq.split('').map((base, index) => {
      const isInRecognitionSite = context_rs_indices.includes(index);
      return (
        <span 
          key={index} 
          className={`font-mono text-lg ${isInRecognitionSite ? 'bg-yellow-200 font-bold' : ''}`}
        >
          {base}
        </span>
      );
    });
  };

  // Render a codon block (original or alternative)
  const renderCodonBlock = (codon, isOriginal = true) => {
    if (!codon) return null;

    // Calculate the position of this codon in the context sequence
    const relativePosition = codon.context_position - context_first_base;
    
    // Style based on whether it's an original or alternative codon
    const blockStyle = {
      position: 'absolute',
      left: `${relativePosition * 10}px`, // Assuming each base is about 10px wide
      top: isOriginal ? '30px' : '0px',
    };

    // Style for bases that overlap with the recognition site
    const getBaseStyle = (baseIndex) => {
      if (isOriginal && codon.rs_overlap && codon.rs_overlap.includes(baseIndex)) {
        return 'bg-yellow-200 font-bold';
      }
      return '';
    };

    return (
      <div style={blockStyle} className="border border-gray-300 px-1 rounded">
        <div className="flex">
          {codon.codon_sequence.split('').map((base, idx) => (
            <span 
              key={idx} 
              className={`font-mono ${getBaseStyle(idx)}`}
            >
              {base}
            </span>
          ))}
        </div>
        {!isOriginal && codon.usage && (
          <div className="text-xs text-center text-gray-600">
            {(codon.usage * 100).toFixed(1)}%
          </div>
        )}
      </div>
    );
  };

  // Function to render the codon comparisons
  const renderCodonComparisons = () => {
    return (
      <div className="mt-12 relative h-48">
        {/* Original codons */}
        {originalCodons.map((codon, index) => (
          <React.Fragment key={`original-${index}`}>
            {renderCodonBlock(codon, true)}
            
            {/* Alternative codon for this position */}
            {renderCodonBlock(
              getAlternativeCodon(index + 1), // +1 because nth_codon_in_rs is 1-based
              false
            )}
          </React.Fragment>
        ))}
      </div>
    );
  };

  // Function to render the mutation set selector
  const renderMutationSetSelector = () => {
    if (!mutationSets || mutationSets.length <= 1) return null;

    return (
      <div className="mb-4">
        <label className="block text-sm font-medium text-gray-700 mb-1">
          Mutation Set:
        </label>
        <select 
          className="border border-gray-300 rounded p-1"
          value={selectedMutationSetIndex}
          onChange={(e) => setSelectedMutationSet(mutationSets[parseInt(e.target.value)])}
        >
          {mutationSets.map((set, idx) => (
            <option key={idx} value={idx}>
              Set {idx + 1} ({set.mutations?.length || 0} mutations)
            </option>
          ))}
        </select>
      </div>
    );
  };

  return (
    <div className="p-4 border rounded shadow-md">
      <h2 className="text-xl font-bold mb-4">Mutation Analysis Summary</h2>
      
      {renderMutationSetSelector()}
      
      <div className="mb-6">
        <h3 className="text-md font-semibold mb-2">Context Sequence with Recognition Site</h3>
        <div className="p-2 border rounded bg-gray-50 overflow-x-auto">
          {renderContextSequence()}
        </div>
        <div className="text-xs mt-1 text-gray-500">
          <span className="inline-block px-1 mx-1 bg-yellow-200">Highlighted</span> bases indicate recognition site
        </div>
      </div>
      
      <div>
        <h3 className="text-md font-semibold mb-2">Codon Replacements</h3>
        <div className="p-2 border rounded bg-gray-50 overflow-x-auto">
          <div className="flex justify-between mb-1">
            <span className="text-xs">Alternative Codons</span>
            <span className="text-xs">Original Codons</span>
          </div>
          {renderCodonComparisons()}
        </div>
        <div className="text-xs mt-1 text-gray-500">
          Percentages indicate codon usage frequency
        </div>
      </div>
    </div>
  );
};

export default MutationAnalysisSummary;