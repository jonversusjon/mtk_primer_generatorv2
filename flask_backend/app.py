import os
import logging
import importlib.util
from flask import Flask, send_from_directory, jsonify, redirect
from flask_cors import CORS

from flask_backend.settings.config import DevelopmentConfig, ProductionConfig
from flask_backend.celery_app import celery_init_app
from flask_backend.routes import api, main

CONFIG_MAP = {
    "Development": DevelopmentConfig,
    "Production": ProductionConfig,
}

def load_prefill_data(env_var="PREFILL_DATA") -> dict:
    path = os.getenv(env_var)
    if path:
        try:
            module_path, attr = path.rsplit(".", 1)
            module = importlib.import_module(module_path)
            return getattr(module, attr)
        except Exception as e:
            logging.error(f"Failed to import dummy data from {path}: {e}")

    # Fallback to default_config
    try:
        from flask_backend.dev.prefill_data.form_defaults import default_config
        logging.info("Using fallback default_config")
        return default_config
    except Exception as e:
        logging.error(f"Failed to load fallback default_config: {e}")
        return {}


def create_app(config_override=None):
    app = Flask(__name__)

    config_path = config_override or os.getenv("CONFIG_OVERRIDE")
    if config_path:
        app.config.from_object(config_path)
    else:
        env = os.getenv("FLASK_ENV", "production")
        app.config.from_object(CONFIG_MAP.get(env.capitalize(), ProductionConfig))
    
    app.config["PREFILL_DATA"] = load_prefill_data()
    
    # Initialize Celery + other extensions
    celery_init_app(app)
    CORS(app, resources={r"/*": {"origins": app.config["CORS_ORIGINS"], "methods": app.config["CORS_METHODS"]}})
    app.register_blueprint(main)
    app.register_blueprint(api, url_prefix="/api")
    
    from flask_backend.routes.sse import custom_sse
    app.register_blueprint(custom_sse, url_prefix="/sse")

    # Redirect /species → /api/species
    @app.route("/species", methods=["GET", "OPTIONS"])
    def species_redirect():
        return redirect("/api/species")

    # Static + React fallbacks
    @app.route("/static/<path:path>")
    def serve_static(path):
        return send_from_directory("static", path)

    @app.route("/", defaults={"path": ""})
    @app.route("/<path:path>")
    def serve_react(path):
        return send_from_directory("static/react", "index.html")

    # JSON error handlers
    @app.errorhandler(404)
    def not_found(error):
        return jsonify(error="Resource not found"), 404

    @app.errorhandler(500)
    def server_error(error):
        return jsonify(error="Internal server error"), 500

    return app
