// src/config/defaultParameters.js
import { templateSequences, testSequences } from "../data/sequences";

// Configuration for different environments/modes
export const appModes = {
  PRODUCTION: "production",
  TESTING: "testing",
  DEVELOPMENT: "development",
};

// Set the current app mode (you can control this with environment variables if needed)
const currentMode = process.env.REACT_APP_MODE || appModes.DEVELOPMENT;

const baseConfig = {};

// Mode-specific configurations
const modeConfigs = {
  [appModes.PRODUCTION]: {
    templateSequence: templateSequences.default,
    sequencesToDomesticate: [],
    mtkPartNums: [],
    primerNames: [],
  },
  [appModes.TESTING]: {
    templateSequence: templateSequences.test,
    sequencesToDomesticate: testSequences,
    mtkPartNums: ["6", "6"],
    primerNames: ["Test 1", "Test 2"],
  },
  [appModes.DEVELOPMENT]: {
    templateSequence: templateSequences.test,
    sequencesToDomesticate: [],
    mtkPartNums: [],
    primerNames: [],
  },
};

export const defaultParameters = {
  ...baseConfig,
  ...modeConfigs[currentMode],
};
