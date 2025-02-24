document.addEventListener("DOMContentLoaded", function () {
    const form = document.getElementById("primer-form");
    if (!form) return;

    const templateSeqEl = document.getElementById("templateSequence");
    
    function isValidDNASequence(sequence) {
        return /^[ATGCWSMKRYBDHVN]*$/i.test(sequence);
    }

    function isAmpliconValid(amplicon, template) {
        if (!template) return true;
        return amplicon.length <= template.length && template.includes(amplicon);
    }

    // HTMX validation handler
    htmx.on("htmx:validation:validate", (evt) => {
        const element = evt.detail.elt;
        const value = element.value.trim();

        // Template sequence validation
        if (element.id === "templateSequence") {
            if (value !== "" && !isValidDNASequence(value)) {
                evt.preventDefault();
                element.setCustomValidity("Only valid DNA bases (A, T, G, C, W, S, M, K, R, Y, B, D, H, V, N) allowed.");
            } else {
                element.setCustomValidity("");
            }
        }

        // Dynamic sequence validation
        if (element.classList.contains("dynamic-sequence-input")) {
            const templateValue = templateSeqEl?.value.trim();

            if (value === "") {
                element.setCustomValidity("Sequence cannot be empty.");
            } else if (!isValidDNASequence(value)) {
                element.setCustomValidity("Only valid DNA bases (A, T, G, C, W, S, M, K, R, Y, B, D, H, V, N) allowed.");
            } else if (value.length < 80) {
                element.setCustomValidity("Amplicon sequence must be at least 80 bp.");
            } else if (templateValue && !isAmpliconValid(value, templateValue)) {
                element.setCustomValidity("Amplicon sequence must be within the template.");
            } else {
                element.setCustomValidity("");
            }

            if (element.validationMessage) {
                evt.preventDefault();
            }
        }

        // MTK Part validation
        if (element.id.startsWith("mtkPart")) {
            if (element.value.trim() === "") {
                evt.preventDefault();
                element.setCustomValidity("MTK Part number is required.");
            } else {
                element.setCustomValidity("");
            }
        }
    });

    // Add validation feedback elements
    htmx.on("htmx:validation:failed", (evt) => {
        const element = evt.detail.elt;
        let feedback = element.parentElement.querySelector('.invalid-feedback');
        
        if (!feedback) {
            feedback = document.createElement('div');
            feedback.className = 'invalid-feedback';
            element.parentElement.appendChild(feedback);
        }
        feedback.textContent = element.validationMessage;
    });

    // Hide validation feedback when field becomes valid
    htmx.on("htmx:validation:validate", (evt) => {
        if (!evt.detail.elt.validationMessage) {
            const feedback = evt.detail.elt.parentElement.querySelector('.invalid-feedback');
            if (feedback) {
                feedback.remove();
            }
        }
    });
});