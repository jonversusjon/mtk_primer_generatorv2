import React from 'react';
import '../../styles/SequenceInput.css';

function SequenceInput({ id, value, onChange, placeholder }) {
  return (
    <div className="sequence-input-wrapper">
      <textarea
        id={id}
        className="sequence-textarea form-control"
        value={value}
        onChange={onChange}
        placeholder={placeholder}
        rows={4}
        spellCheck="false"
      />
    </div>
  );
}

export default SequenceInput;