from flask_cors import CORS
from config.logging_config import logger
from flask import Flask, send_from_directory
from routes.main import main
from routes.api import api
from config.settings import Config, TestConfig
import os


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
    app.config.from_object(config_class)
    logger.info(f"TESTING mode: {app.config['TESTING']}")

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
        # You can either redirect or implement the same functionality here
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
        return send_from_directory('static', path)

    # Serve React app - in production, this would be handled by a web server
    @app.route('/', defaults={'path': ''})
    @app.route('/<path:path>')
    def serve_react(path):
        return send_from_directory('static/react', 'index.html')

    # Error handlers
    @app.errorhandler(404)
    def not_found(error):
        return {'error': 'Resource not found'}, 404

    @app.errorhandler(500)
    def server_error(error):
        return {'error': 'Internal server error'}, 500

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(debug=True)
