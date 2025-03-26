export const updateFromSSE = (prevState, sseData) => ({
    ...prevState,
    ...sseData,
  });
  