from flask_cors import CORS
from config.logging_config import logger
from flask import Flask, send_from_directory
from flask import Flask
from routes.main import main
from routes.api import api
from routes.validation import validation
from config.settings import Config, TestConfig
import os
from config.logging_config import logger  # Use centralized logger


def configure_werkzeug_logging():
    """Configure Werkzeug to use our logging settings."""
    werkzeug_logger = logger.getChild("werkzeug")
    werkzeug_logger.setLevel(logger.level)


def create_app():
    """Create and configure the Flask app."""
    testing_env = os.getenv('FLASK_TESTING', 'false').lower() == 'true'
    print(f"FLASK_TESTING environment variable: {testing_env}")

    config_class = TestConfig if testing_env else Config

    app = Flask(__name__)
    app.config.from_object(config_class)
    print(f"TESTING mode: {app.config['TESTING']}")

    logger.info(f"Starting Flask app with config: {config_class.__name__}")

    configure_werkzeug_logging()

    app.register_blueprint(main)
    app.register_blueprint(validation, url_prefix='/validation')
    app.register_blueprint(api, url_prefix="/api")

    @app.context_processor
    def utility_processor():
        """Make test settings available to Jinja templates."""
        return {
            'TESTING_MODE': app.config['TESTING'],
            'TEST_SEQ': app.config.get('TEST_SEQ', ''),
            'TEST_TEMPLATE_SEQ': app.config.get('TEST_TEMPLATE_SEQ', '')
        }

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(debug=True)


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
    CORS(app)

    logger.info(f"Starting Flask app with config: {config_class.__name__}")

    configure_werkzeug_logging()

    # Register API blueprint
    app.register_blueprint(api)

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
