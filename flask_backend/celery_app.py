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
    celery.conf.update(app.config["CELERY"])
    celery.Task = FlaskTask
    # Save the celery instance in app.extensions for later access.
    app.extensions["celery"] = celery
    return celery
