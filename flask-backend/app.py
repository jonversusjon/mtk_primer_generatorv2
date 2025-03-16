import os
import importlib.util

import argparse
from flask import Flask, send_from_directory, jsonify
from flask_cors import CORS
from config.logging_config import logger
from routes.main import main
from routes.api import api
from config.settings import Config, TestConfig


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
            # Print full config
            # print(f"🔍 DEBUG: Full CONFIG dictionary: {config_data}")

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


def configure_werkzeug_logging():
    """Configure Werkzeug to use our logging settings."""
    werkzeug_logger = logger.getChild("werkzeug")
    werkzeug_logger.setLevel(logger.level)


def create_app():
    """Create and configure the Flask app."""
    testing_env = os.getenv('FLASK_TESTING', 'false').lower() == 'true'
    logger.info(f"FLASK_TESTING environment variable: {testing_env}")

    config_class = TestConfig if testing_env else Config

    app = Flask(__name__)
    app.config["ACTIVE_CONFIG"] = app_config
    # logger.info(f"TESTING mode: {app.config['TESTING']}")
    # logger.info(f"🔍 DEBUG: Loaded config: {app.config['ACTIVE_CONFIG']}")

    # Enable CORS for the entire app
    CORS(app, resources={
        r"/*": {
            "origins": ["http://localhost:3000"],  # React dev server
            "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"],
            "allow_headers": ["Content-Type"]
        }
    })
    logger.info(f"Starting Flask app with config: {config_class.__name__}")

    configure_werkzeug_logging()

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

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(debug=True)
