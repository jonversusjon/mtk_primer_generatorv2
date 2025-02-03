from flask import Blueprint, render_template, request, jsonify
from tests.test_data import TEST_SEQ, TEST_TEMPLATE_SEQ
from utils.utils import get_available_species
from services.protocol import create_gg_protocol

main_blueprint = Blueprint("main", __name__)

TESTING_MODE = True

@main_blueprint.route("/")
def home():

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

@main_blueprint.route('/validate_inputs', methods=['POST'])
def validate_inputs():
    # Your request handling logic here...
    # Get data from request.get_json()
    # Process the data...
    # Return a response (e.g., jsonify({'message': 'Success'}))
    return "Placeholder"

@main_blueprint.route('/generate_protocol', methods=['POST'])
def generate_protocol():
    """
    Generates the protocol for primer design based on the number of sequences provided in the request.

    The function parses the input form data, which includes the number of sequences, species, Kozak sequence type,
    and verbose mode status. It then iterates over each sequence, processes it, and generates the corresponding
    primer design protocol.

    Returns:
        str: A string containing the generated protocol for primer design or an error message.
    """

    if request.method == 'POST':
        try:
            data = request.form.to_dict(flat=False)  # Parse form data to dictionary
            print(f'data: {data}')
            # Extract the number of sequences and initialize the response
            num_sequences = int(data.get('numSequences', [0])[0])
            if num_sequences == 0:
                return "No sequences provided.", 400

            kozak = data.get('kozak', [''])[0]
            species = data.get('species', [''])[0]
            verbose = 'verbose_mode' in data

            print(f'kozak: {kozak}, species: {species}, verbose: {verbose}')

            response = ""
            for i in range(num_sequences):
                seq_data = data.get(f'sequences[{i}][]', [])
                if len(seq_data) < 3:
                    return f"Incomplete data for sequence {i+1}.", 400
                
                primer_name = seq_data[0]
                mtk_part = seq_data[1]
                sequence = seq_data[2]

            #     # Process each sequence
            #     protocol = process_sequence(
            #         primer_name,
            #         mtk_part,
            #         sequence,
            #         species,
            #         kozak,
            #         verbose=verbose
            #     )
                print(f"Processing sequence {i+1}...\n")
                # response += f"Protocol for sequence {i+1}:\n{protocol}\n\n"

            return response, 200
        except Exception as e:
            print(f"Error processing request: {e}")
            return "An error occurred while processing the request.", 500

    return "Invalid request method.", 405