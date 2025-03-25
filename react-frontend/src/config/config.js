// Base URL for API requests - adjust based on your environment
export const SSE_BASE_URL = 
  (process.env.REACT_APP_SSE_URL || "http://localhost:5000") + "/sse";

export const API_BASE_URL =
  (process.env.REACT_APP_API_URL || "http://localhost:5000") + "/api";
