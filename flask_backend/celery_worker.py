import flask_backend.celery_tasks # noqa

from .app import create_app

# Create the Flask app and initialize Celery
flask_app = create_app()
celery = flask_app.extensions["celery"]
