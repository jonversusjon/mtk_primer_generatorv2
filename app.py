from flask import Flask
from routes.main import main
from routes.api import api_blueprint
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
    app.register_blueprint(api_blueprint, url_prefix="/api")

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
