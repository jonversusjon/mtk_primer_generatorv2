from flask import Flask
from routes.main import main
from routes.api import api_blueprint
from config.settings import Config, TestConfig
import os


def create_app():
    testing_env = os.getenv('FLASK_TESTING', '').lower()
    print(f"FLASK_TESTING environment variable: {testing_env}")

    config_class = TestConfig if testing_env == 'true' else Config
    print(f"Using config class: {config_class.__name__}")

    app = Flask(__name__)
    app.config.from_object(config_class)
    
    # Print relevant config values
    print(f"TESTING mode: {app.config['TESTING']}")
    print(f"TEST_SEQ available: {'TEST_SEQ' in app.config}")
    print(f"TEST_TEMPLATE_SEQ available: {'TEST_TEMPLATE_SEQ' in app.config}")

    @app.context_processor
    def utility_processor():
        return {
            'TESTING_MODE': app.config['TESTING'],
            'TEST_SEQ': app.config.get('TEST_SEQ', ''),
            'TEST_TEMPLATE_SEQ': app.config.get('TEST_TEMPLATE_SEQ', '')
        }

    app.register_blueprint(main)
    app.register_blueprint(api_blueprint, url_prefix="/api")

    return app

