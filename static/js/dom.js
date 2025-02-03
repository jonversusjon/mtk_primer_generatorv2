// dom.js

function updateSequenceInputs(num, defaultTestSeq = "") {
    let container = document.getElementById("sequence-container");
    num = Math.max(1, Math.min(10, num)); // Ensure number stays within range
    container.innerHTML = "";

    // Create tabs for the sequences
    let tabs = document.createElement("div");
    tabs.classList.add("sequence-tabs");
    for (let i = 0; i < num; i++) {
        let tab = document.createElement("button");
        tab.classList.add("sequence-tab");
        tab.textContent = `Sequence ${i + 1}`;
        tab.addEventListener("click", function () {
            openSequenceTab(i);
        });
        if (i === 0) {
            tab.classList.add("active");
        }
        tabs.appendChild(tab);
    }
    container.appendChild(tabs);

    // MTK Part Number options
    const mtk_part_nums = [
        '', '1', '2', '3', '3a', '3b', '3c', '3d', '3e', '3f',
        '3g', '4', '4a', '4b', '4aII', '3bI', '3bII', '5', '6', '7', '8', '8a', '8b'
    ];

    // Create content panes for each sequence
    let content = document.createElement("div");
    content.classList.add("sequence-content");
    for (let i = 0; i < num; i++) {
        let pane = document.createElement("div");
        pane.classList.add("sequence-pane");
        if (i === 0) {
            pane.classList.add("active");
        }

        // Prepopulate test values only if `TESTING_MODE` is true
        var TESTING_MODE = JSON.parse(document.getElementById("testing-mode-data").textContent);

        let testPrimerName = TESTING_MODE && i === 0 ? "Test Primer 1" : "";
        let testSequence = TESTING_MODE && i === 0 ? defaultTestSeq : "";
        let testPartNumber = TESTING_MODE && i === 0 ? "6" : "";

        pane.innerHTML = `
            <div class="form-row-top">
                <div class="input-group">
                    <label for="primerName${i}">Primer Name:</label>
                    <input type="text" id="primerName${i}" name="sequences[${i}][primerName]" class="form-control"
                        placeholder="Enter primer name (optional)"
                        value="${testPrimerName}" />
                </div>

                <div class="input-group">
                    <label for="mtkPart${i}">MTK Part Number:</label>
                    <select id="mtkPart${i}" name="sequences[${i}][mtkPart]" class="form-select">
                        ${mtk_part_nums.map(val =>
            `<option value="${val}" ${val === testPartNumber ? "selected" : ""}>${val}</option>`
        ).join('')}
                    </select>
                </div>
            </div>

            <label>Sequence ${i + 1}:</label>
            <textarea id="sequence${i}" name="sequences[${i}][sequence]" class="form-control dynamic-sequence-input sequence-input"
                placeholder="Enter DNA sequence...">${testSequence}</textarea>
            <div class="char-count-label" style="display:none;"></div>
        `;

        content.appendChild(pane);

        // Attach event listener to dynamically created sequence inputs for char count
        const sequenceInput = pane.querySelector('.dynamic-sequence-input');
        sequenceInput.addEventListener('input', function () {
            updateCharCount(sequenceInput);
        });

        // Trigger 'input' event if defaultTestSeq is provided
        if (testSequence) {
            sequenceInput.dispatchEvent(new Event('input'));
        }
    }
    container.appendChild(content);
}

document.addEventListener('DOMContentLoaded', function () {
    // Get the template sequence element
    const templateSeqEl = document.getElementById('templateSequence');

    // Function to update character count for a given element
    function updateCharCount(el) {
        let label = el.nextElementSibling;
        if (!label || !label.classList.contains("char-count-label")) return;
        let len = el.value.length;
        if (len > 0) {
            label.style.display = "block";
            label.textContent = `Length: ${len} bp`;
        } else {
            label.style.display = "none";
        }
    }

    // Update char count for template sequence if it has a default value
    if (templateSeqEl && templateSeqEl.value) {
        updateCharCount(templateSeqEl);
    }

    // Attach event listener to template sequence input for future changes
    if (templateSeqEl) {
        templateSeqEl.addEventListener('input', function () {
            updateCharCount(templateSeqEl);
        });
    }

    // Initialize dynamic sequence inputs
    updateSequenceInputs(document.getElementById("numSequences").value, defaultTestSeq);

    // Update char count for dynamically created sequence inputs if they have default values
    document.querySelectorAll('.dynamic-sequence-input').forEach(input => {
        if (input.value) {
            updateCharCount(input);
        }
    });
});

function updateCharCount(el) {
    let label = el.nextElementSibling;
    if (!label || !label.classList.contains("char-count-label")) return;
    let len = el.value.length;
    if (len > 0) {
        label.style.display = "block";
        label.textContent = `Length: ${len} bp`;
    } else {
        label.style.display = "none";
    }
}

// Function to switch between sequence panes
function openSequenceTab(index) {
    document.querySelectorAll('.sequence-tab').forEach(tab => tab.classList.remove('active'));
    document.querySelectorAll('.sequence-pane').forEach(pane => pane.classList.remove('active'));
    document.querySelectorAll('.sequence-tab')[index].classList.add('active');
    document.querySelectorAll('.sequence-pane')[index].classList.add('active');
}

// Clear the form content
function clearForm() {
    document.getElementById("fileUpload").value = "";
    document.getElementById("numSequences").value = 1;
    document.getElementById("sequence-container").innerHTML = "";
    document.getElementById("results").innerHTML = "";
}

// Tab switching functionality for Template Sequence
function openTab(tabName) {
    const tabPanes = document.querySelectorAll('.tab-pane');
    tabPanes.forEach(pane => pane.classList.remove('active'));
    const tabButtons = document.querySelectorAll('.tab-button');
    tabButtons.forEach(button => button.classList.remove('active'));
    document.getElementById(tabName).classList.add('active');
    event.currentTarget.classList.add('active');
}