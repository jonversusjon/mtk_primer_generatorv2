import os
import json
from flask import Blueprint, render_template, request
from tests.test_data import TEST_SEQ, TEST_TEMPLATE_SEQ
from utils.utils import get_available_species

main_blueprint = Blueprint("main", __name__)

TESTING_MODE = True

@main_blueprint.route("/")
def home():
    print(f"DEBUG: test_seq = {repr(TEST_SEQ)[:100]}...")  # Print first 100 chars
    print(f"DEBUG: test_template_seq = {repr(TEST_TEMPLATE_SEQ)[:100]}...")
    return render_template(
        "index.html",
        title="Home Page",
        test_seq=TEST_SEQ,
        test_template_seq=TEST_TEMPLATE_SEQ,
        testing_mode=TESTING_MODE
    )

@main_blueprint.route("/get_species")
def get_species():
    species = get_available_species()
    options_html = "".join(f'<option value="{s}">{s}</option>' for s in species)
    return options_html  # Return only the options, not the <select>


@main_blueprint.route("/get_sequence_inputs", methods=["GET"])
def get_sequence_inputs():
    num_sequences = int(request.args.get("numSequences", 1))  # Get the user-inputted number of sequences
    inputs_html = ""

    for i in range(num_sequences):
        inputs_html += f"""
        <div class="sequence-input-group">
            <label for="sequence-{i}" class="form-label">Sequence {i+1}:</label>
            <textarea id="sequence-{i}" name="sequences[]" class="form-control sequence-input"
                      placeholder="Enter DNA sequence here..."></textarea>
            <div id="charCount-sequence-{i}" class="char-count-label" style="display:none;"></div>
        </div>
        """

    return inputs_html  # Return the dynamically generated input fields
