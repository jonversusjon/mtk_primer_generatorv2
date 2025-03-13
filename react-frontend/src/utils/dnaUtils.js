/**
 * Utility functions for DNA sequence validation and manipulation
 */

/**
 * Validates a DNA sequence
 * @param {string} sequence - The DNA sequence to validate
 * @param {boolean} isRequired - Whether the sequence is required
 * @param {boolean} requireInFrame - Whether the sequence must be divisible by 3
 * @returns {Object} Object with isValid status, message, and isAdvisory flag
 */
export const validateDnaSequence = (
  sequence,
  isRequired = false,
  requireInFrame = true
) => {
  // Check for empty sequence if required
  const trimmedSequence = typeof sequence === "string" ? sequence.trim() : "";

  if (isRequired && !trimmedSequence) {
    return {
      isValid: false,
      message: "Sequence cannot be empty",
      isAdvisory: false,
    };
  }

  if (!sequence || sequence.trim() === "") {
    return {
      isValid: true,
      message: "",
      isAdvisory: false,
    };
  }

  // Clean the sequence (remove whitespace)
  const cleanedSequence = sequence.trim().toUpperCase();

  // Check for valid DNA bases (A, T, G, C, plus extended DNA codes)
  const validDnaRegex = /^[ATGCWSMKRYBDHVN]+$/;
  if (!validDnaRegex.test(cleanedSequence)) {
    return {
      isValid: false,
      message:
        "Only valid DNA bases (A, T, G, C) or extended DNA code (W, S, M, K, R, Y, B, D, H, V, N) allowed",
      isAdvisory: false,
    };
  }

  // Check length (minimum 80 bp for meaningful sequences)
  if (cleanedSequence.length < 80) {
    return {
      isValid: false,
      message: "Sequence must be at least 80 bp",
      isAdvisory: false,
    };
  }

  // Check if sequence is in frame (divisible by 3) - only if required
  if (requireInFrame) {
    const inFrame = cleanedSequence.length % 3 === 0;
    if (!inFrame) {
      return {
        isValid: false,
        message:
          "Sequence appears to be out of frame (not divisible by 3). Double-check your sequence if this is meant to be entirely a coding sequence.",
        isAdvisory: true,
      };
    }
  }

  // Check for start and stop codons (only if we're requiring in-frame)
  if (requireInFrame) {
    const hasStartCodon = cleanedSequence.substring(0, 3) === "ATG";
    const stopCodons = ["TAA", "TAG", "TGA"];
    const lastCodon = cleanedSequence.substring(cleanedSequence.length - 3);
    const hasStopCodon = stopCodons.includes(lastCodon);

    // Create advisory messages for codons
    if (hasStartCodon && hasStopCodon) {
      return {
        isValid: true,
        message:
          "Start and stop codons detected and will be removed during processing",
        isAdvisory: true,
      };
    } else if (hasStartCodon) {
      return {
        isValid: true,
        message: "Start codon detected and will be removed during processing",
        isAdvisory: true,
      };
    } else if (hasStopCodon) {
      return {
        isValid: true,
        message: "Stop codon detected and will be removed during processing",
        isAdvisory: true,
      };
    }
  }

  // Sequence is valid with no advisories
  return {
    isValid: true,
    message: "",
    isAdvisory: false,
  };
};

/**
 * Reverse complements a DNA sequence
 * @param {string} sequence - The DNA sequence to reverse complement
 * @returns {string} The reverse complemented sequence
 */
export const reverseComplement = (sequence) => {
  const complement = {
    A: "T",
    T: "A",
    G: "C",
    C: "G",
    a: "t",
    t: "a",
    g: "c",
    c: "g",
    W: "W",
    S: "S",
    M: "K",
    K: "M",
    R: "Y",
    Y: "R",
    B: "V",
    D: "H",
    H: "D",
    V: "B",
    N: "N",
    w: "w",
    s: "s",
    m: "k",
    k: "m",
    r: "y",
    y: "r",
    b: "v",
    d: "h",
    h: "d",
    v: "b",
    n: "n",
  };

  return sequence
    .split("")
    .reverse()
    .map((base) => complement[base] || base)
    .join("");
};

/**
 * Calculates GC content percentage of a DNA sequence
 * @param {string} sequence - The DNA sequence
 * @returns {number} GC content as a percentage
 */
export const calculateGCContent = (sequence) => {
  if (!sequence) return 0;

  const cleanedSequence = sequence.toUpperCase();
  const gcCount = (cleanedSequence.match(/[GC]/g) || []).length;

  return (gcCount / cleanedSequence.length) * 100;
};

/**
 * Calculates the approximate melting temperature of a DNA sequence
 * @param {string} sequence - The DNA sequence
 * @returns {number} Approximate melting temperature in °C
 */
export const calculateMeltingTemp = (sequence) => {
  if (!sequence) return 0;

  const cleanedSequence = sequence.toUpperCase();
  const length = cleanedSequence.length;

  // Count nucleotides
  const aCount = (cleanedSequence.match(/A/g) || []).length;
  const tCount = (cleanedSequence.match(/T/g) || []).length;
  const gCount = (cleanedSequence.match(/G/g) || []).length;
  const cCount = (cleanedSequence.match(/C/g) || []).length;

  // Basic melting temperature calculation
  let meltingTemp = 0;

  if (length < 14) {
    // Formula for short sequences
    meltingTemp = (aCount + tCount) * 2 + (gCount + cCount) * 4;
  } else {
    // Wallace formula for longer sequences
    meltingTemp = 64.9 + (41 * (gCount + cCount - 16.4)) / length;
  }

  return Math.round(meltingTemp * 100) / 100; // Round to 2 decimal places
};
