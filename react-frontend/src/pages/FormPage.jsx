import React, { useState, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import Form from "../components/Form/Form";
import { generateProtocol } from "../api/api";
import "../styles/Form.css";

function FormPage({ showSettings, setShowSettings, setResults }) {
  const [formData, setFormData] = useState({ sequencesToDomesticate: [""] });
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(null);

  const navigate = useNavigate();

  useEffect(() => {
    const fetchConfig = async () => {
      try {
        const savedData = sessionStorage.getItem("formData");

        if (savedData) {
          setFormData(JSON.parse(savedData));
        } else {
          const response = await fetch(`http://localhost:5000/api/config`);
          const data = await response.json();

          if (!data || Object.keys(data).length === 0) {
            console.warn("API returned empty config! Using fallback defaults.");
            setFormData({ sequencesToDomesticate: [""] });
          } else {
            setFormData({
              ...data,
              sequencesToDomesticate: data.sequencesToDomesticate || [""], // Ensure it's always an array
            });
          }
        }
      } catch (err) {
        console.error("Error fetching defaults from Flask:", err);
        setFormData({ sequencesToDomesticate: [""] }); // Fallback in case of error
      } finally {
        setLoading(false);
      }
    };

    fetchConfig();
  }, []);

  const handleFormSubmit = async (data) => {
    setLoading(true);
    setError(null);

    try {
      sessionStorage.setItem("formData", JSON.stringify(data));
      const response = await generateProtocol(data);
      sessionStorage.setItem("results", JSON.stringify(response));
      setResults(response);
      navigate("/results");
    } catch (err) {
      setError(
        err.message || "An error occurred while generating the protocol"
      );
      console.error("Error generating protocol:", err);
    } finally {
      setLoading(false);
    }
  };

  if (loading) {
    return <p>Loading defaults...</p>;
  }

  return (
    <div className="form-container">
      {error && <div className="alert alert-danger">{error}</div>}

      <Form
        onSubmit={handleFormSubmit}
        onChange={setFormData}
        formData={formData ?? { sequencesToDomesticate: [""] }}
        showSettings={showSettings}
        setShowSettings={setShowSettings}
      />
    </div>
  );
}

export default FormPage;
