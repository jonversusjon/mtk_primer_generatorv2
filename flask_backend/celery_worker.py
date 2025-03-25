from flask_backend.app import create_app

# Build the Flask app (this reads your ENV vars, including CELERY_BROKER_URL)
flask_app = create_app()

# Grab the Celery instance that create_app() already set up.
celery = flask_app.extensions["celery"]

# Import tasks so that Celery can auto-discover them.
import flask_backend.celery_tasks  # noqa
