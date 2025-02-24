from flask import Blueprint, jsonify
from tests.test_data import TEST_SEQ, TEST_TEMPLATE_SEQ

api_blueprint = Blueprint("api", __name__)

# @api_blueprint.route("/data")
# def get_data():
#     return jsonify({"message": "Hello, API!"})

# @api_blueprint.route("/test-sequence")
# def get_test_seq():
#     return {"test_sequence": TEST_SEQ}
            