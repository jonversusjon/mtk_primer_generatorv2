import os
import importlib.util
from flask import Flask, send_from_directory, jsonify, request, stream_with_context, Response
from celery_worker import celery
import time
import json
import argparse
from flask_cors import CORS
from routes.main import main
from routes.api import api
from config.settings import Config, TestConfig
from log_utils import logger

def load_python_config(module_path, env="development"):
    """Dynamically load a Python config module and return the correct environment settings."""
    try:
        spec = importlib.util.find_spec(module_path)
        if spec is None:
            raise ImportError(f"Module '{module_path}' not found.")

        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        if hasattr(module, "CONFIG"):
            config_data = module.CONFIG

            # If CONFIG contains multiple environments, extract the correct one
            if isinstance(config_data, dict) and env in config_data:
                selected_config = config_data[env]
                # print(f"✅ DEBUG: Extracted '{env}' config: {selected_config}")
                return selected_config
            return config_data  # If it's not environment-based, return as is
        else:
            raise AttributeError(
                f"Module '{module_path}' does not contain a 'CONFIG' dictionary.")
    except Exception as e:
        logger.error(f"Error loading config module '{module_path}': {e}")
        return {}


# Argument parsing for config file
parser = argparse.ArgumentParser(
    description="Start the Flask app with a custom config file.")
parser.add_argument("--config", type=str, default="config.default_config",
                    help="Path to the config module (dot notation).")
parser.add_argument("--env", type=str, default="development",
                    help="Configuration environment (development/testing/production).")
args = parser.parse_args()

CONFIG_MODULE = args.config
ENVIRONMENT = args.env

# print(f"🔍 DEBUG: Loading {CONFIG_MODULE} with environment '{ENVIRONMENT}'")

app_config = load_python_config(CONFIG_MODULE, ENVIRONMENT)


def create_app():
    """Create and configure the Flask app."""
    testing_env = os.getenv('FLASK_TESTING', 'false').lower() == 'true'
    logger.log_step("", f"FLASK_TESTING environment variable: {testing_env}")

    config_class = TestConfig if testing_env else Config

    app = Flask(__name__)
    app.config["ACTIVE_CONFIG"] = app_config
    
    # Enable CORS for the entire app
    CORS(app, resources={
        r"/*": {
            "origins": ["http://localhost:3000"],  # React dev server
            "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"],
            "allow_headers": ["Content-Type"]
        }
    })
    logger.log_step("", f"Starting Flask app with config: {config_class.__name__}")

    # Register blueprints
    app.register_blueprint(main)
    app.register_blueprint(api, url_prefix="/api")

    @app.route('/species', methods=['GET', 'OPTIONS'])
    def species_redirect():
        """Redirect /species to /api/species."""
        from flask import redirect
        return redirect('/api/species')

    # Make test settings available to Jinja templates
    @app.context_processor
    def utility_processor():
        return {
            'TESTING_MODE': app.config['TESTING'],
            'TEST_SEQ': app.config.get('TEST_SEQ', ''),
            'TEST_TEMPLATE_SEQ': app.config.get('TEST_TEMPLATE_SEQ', '')
        }

    # Serve static files
    @app.route('/static/<path:path>')
    def serve_static(path):
        """Serve static files."""
        return send_from_directory('static', path)

    # Serve React app - in production, this would be handled by a web server
    @app.route('/', defaults={'path': ''})
    @app.route('/<path:path>')
    def serve_react(path):
        """Serve React frontend files."""
        return send_from_directory('static/react', 'index.html')

    # Error handlers
    @app.errorhandler(404)
    def not_found(error):
        """Handle 404 errors with JSON response."""
        return jsonify({'error': 'Resource not found'}), 404

    @app.errorhandler(500)
    def server_error(error):
        """Handle 500 errors with JSON response."""
        return jsonify({'error': 'Internal server error'}), 500

    @app.route("/api/start-task", methods=["POST"])
    def start_task():
        """Enqueue a Celery task; payload determines which function to call."""
        payload = request.get_json()
        task = celery.send_task(payload["task_name"], args=payload.get("args", []), kwargs=payload.get("kwargs", {}))
        return jsonify({"task_id": task.id}), 202

    @app.route("/api/task-status/<task_id>")
    def task_status(task_id):
        def event_stream():
            while True:
                res = celery.AsyncResult(task_id)
                status = res.status
                data = {"task_id": task_id, "status": status}
                if res.status == "SUCCESS":
                    data["result"] = res.result
                    yield f"data: {json.dumps(data)}\n\n"
                    break
                elif res.status == "FAILURE":
                    data["error"] = str(res.result)
                    yield f"data: {json.dumps(data)}\n\n"
                    break
                else:
                    yield f"data: {json.dumps(data)}\n\n"
                time.sleep(1)
        return Response(stream_with_context(event_stream()), mimetype="text/event-stream")

    return app

if __name__ == "__main__":
    app = create_app()
    app.run(debug=True, host="0.0.0.0", port=5000)
