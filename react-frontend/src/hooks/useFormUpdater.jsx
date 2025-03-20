import { useCallback } from "react";

export const useFormUpdater = (setFormData) => {
  // For updates coming from the Settings panel
  const updateSettings = useCallback(
    (keyOrUpdate, value) => {
      setFormData((prev) => {
        if (typeof keyOrUpdate === "string" && value !== undefined) {
          return { ...prev, [keyOrUpdate]: value };
        } else {
          const updateObj =
            typeof keyOrUpdate === "function" ? keyOrUpdate(prev) : keyOrUpdate;
          return { ...prev, ...updateObj };
        }
      });
    },
    [setFormData]
  );

  // For updates coming from the Form inputs
  const updateFormInput = useCallback(
    (keyOrUpdate, value) => {
      setFormData((prev) => {
        if (typeof keyOrUpdate === "string" && value !== undefined) {
          return { ...prev, [keyOrUpdate]: value };
        } else {
          const updateObj =
            typeof keyOrUpdate === "function" ? keyOrUpdate(prev) : keyOrUpdate;
          return { ...prev, ...updateObj };
        }
      });
    },
    [setFormData]
  );

  return { updateSettings, updateFormInput };
};
