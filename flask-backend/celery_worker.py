from celery import Celery
from app import create_app

flask_app = create_app()
celery = Celery(
    flask_app.import_name,
    broker=flask_app.config.get("CELERY_BROKER_URL"),
    backend=flask_app.config.get("CELERY_RESULT_BACKEND")
)
celery.conf.update(flask_app.config)

# Ensure tasks run with Flask context
TaskBase = celery.Task
class ContextTask(TaskBase):
    def __call__(self, *args, **kwargs):
        with flask_app.app_context():
            return TaskBase.__call__(self, *args, **kwargs)

celery.Task = ContextTask
