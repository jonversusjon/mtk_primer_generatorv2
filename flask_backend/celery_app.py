import logging
from celery import Celery, Task

def celery_init_app(app):
    """
    Initialize Celery with Flask’s configuration.
    This sets up a FlaskTask so that all tasks run within the app context.
    """
    class FlaskTask(Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)
    celery = Celery(app.import_name)
    
    # Disable Celery logging
    logging.getLogger('celery').propagate = False
    logging.getLogger('celery').handlers = []

    celery.conf.update(app.config["CELERY"])
    
    celery.conf.update(
        worker_hijack_root_logger=False,
        worker_redirect_stdouts=False
    )
    
    celery.Task = FlaskTask
    # Save the celery instance in app.extensions for later access.
    app.extensions["celery"] = celery
    return celery
