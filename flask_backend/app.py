import os
import importlib.util
from flask import Flask, send_from_directory, jsonify, redirect
from flask_cors import CORS

from flask_backend.celery_app import celery_init_app
from flask_backend.routes import api, main
from flask_backend.log_utils import logger

def load_python_config(module_path):
    try:
        spec = importlib.util.find_spec(module_path)
        if spec is None:
            raise ImportError(f"Module '{module_path}' not found.")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        if hasattr(module, "CONFIG"):
            return module.CONFIG
        else:
            raise AttributeError(f"Module '{module_path}' does not contain a 'CONFIG' dictionary.")
    except Exception as e:
        logger.error(f"Error loading config module '{module_path}': {e}")
        return {}


def create_app(config_module=None, env="development"):
    """Create and configure the Flask app."""

    if not config_module:
        config_module = os.getenv("CONFIG_MODULE")
        
    app_config = load_python_config(config_module) if config_module else {}

    app = Flask(__name__)
    celery_init_app(app)

    # Enable CORS for the entire app
    CORS(app, resources={
        r"/*": {
            "origins": ["http://localhost:3000"],
            "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"],
            "allow_headers": ["Content-Type"]
        }
    })
    
    # Register blueprints
    app.register_blueprint(main)
    app.register_blueprint(api, url_prefix="/api")
    
    app.config["ACTIVE_CONFIG"] = app_config
    print(f"Loaded config: {app_config}")
    
    # Inject CELERY configuration into Flask config.
    app.config["CELERY"] = {
        "broker_url": os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/0"),
        "result_backend": os.getenv("CELERY_RESULT_BACKEND", "redis://localhost:6379/0"),
        "task_serializer": "json",
        "result_serializer": "json",
        "accept_content": ["json"],
        "task_ignore_result": False,
    }

    @app.route('/species', methods=['GET', 'OPTIONS'])
    def species_redirect():
        """Redirect /species to /api/species."""
        return redirect('/api/species')

    # Make test settings available to Jinja templates
    @app.context_processor
    def utility_processor():
        return {
            'TESTING_MODE': app.config.get('TESTING', False),
            'TEST_SEQ': app.config.get('TEST_SEQ', ''),
            'TEST_TEMPLATE_SEQ': app.config.get('TEST_TEMPLATE_SEQ', '')
        }

    # Serve static files
    @app.route('/static/<path:path>')
    def serve_static(path):
        """Serve static files."""
        return send_from_directory('static', path)

    # Serve React app (in production, typically handled by a web server)
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
    import argparse

    parser = argparse.ArgumentParser(
        description="Start the Flask app with a custom config file."
    )
    parser.add_argument("--config", type=str, default="config.default_config",
                        help="Path to the config module (dot notation).")
    parser.add_argument("--env", type=str, default="development",
                        help="Configuration environment (development/testing/production).")
    args = parser.parse_args()

    app = create_app(config_module=args.config, env=args.env)
    app.run(debug=True, host="0.0.0.0", port=5000)
