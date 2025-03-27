export function isObject(item) {
  return item !== null && typeof item === "object" && !Array.isArray(item);
}

function isEmpty(value) {
  if (value === null || value === undefined) return true;
  if (typeof value === "string" && value.trim() === "") return true;
  if (Array.isArray(value) && value.length === 0) return true;
  if (isObject(value) && Object.keys(value).length === 0) return true;
  return false;
}

/**
 * Recursively merges source into target.
 * - Only overwrites target if it is null/empty.
 * - If target already has a non-empty value:
 *    • For arrays: concatenates source array.
 *    • For objects: recursively merges.
 *    • For strings: appends source string (if not already present).
 *    • Otherwise (numbers, booleans, mismatched types): leaves target unchanged.
 */
export function mergeDeep(target, source) {
  if (!source) return target;
  const output = { ...target };

  for (const key in source) {
    const sourceVal = source[key];
    const targetVal = output[key];

    // If target is empty, use the source value.
    if (isEmpty(targetVal)) {
      output[key] = sourceVal;
      continue;
    }

    // If both are arrays, concatenate them.
    if (Array.isArray(targetVal) && Array.isArray(sourceVal)) {
      output[key] = targetVal.concat(sourceVal);
      continue;
    }

    // If both are objects, merge recursively.
    if (isObject(targetVal) && isObject(sourceVal)) {
      output[key] = mergeDeep(targetVal, sourceVal);
      continue;
    }

    // If both are strings, append new data if not already contained.
    if (typeof targetVal === "string" && typeof sourceVal === "string") {
      if (!targetVal.includes(sourceVal)) {
        output[key] = targetVal + " " + sourceVal;
      }
      continue;
    }

    // For numbers, booleans, or mismatched types, leave the target value unchanged.
  }

  return output;
}
