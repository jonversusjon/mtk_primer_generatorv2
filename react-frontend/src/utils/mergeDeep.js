// mergeDeep.js
export function isObject(item) {
  return item !== null && typeof item === "object" && !Array.isArray(item);
}

/**
 * Recursively merges source into target.
 * - For objects, it performs a deep merge.
 * - For arrays, it concatenates them.
 * - For primitives, the source overwrites target.
 *
 * @param {object} target - The original object.
 * @param {object} source - The object with new values.
 * @returns {object} A new object with merged properties.
 */
export function mergeDeep(target, source) {
  if (!source) return target;
  const output = { ...target };

  for (const key in source) {
    const sourceVal = source[key];
    const targetVal = output[key];

    if (Array.isArray(sourceVal)) {
      // Concatenate arrays (you can change this behavior as needed)
      if (Array.isArray(targetVal)) {
        output[key] = [...targetVal, ...sourceVal];
      } else {
        output[key] = [...sourceVal];
      }
    } else if (isObject(sourceVal)) {
      if (!targetVal) {
        output[key] = { ...sourceVal };
      } else if (isObject(targetVal)) {
        output[key] = mergeDeep(targetVal, sourceVal);
      } else {
        output[key] = { ...sourceVal };
      }
    } else {
      output[key] = sourceVal;
    }
  }

  return output;
}
