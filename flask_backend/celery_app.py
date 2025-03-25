from celery import Celery, Task

def celery_init_app(app):
    class FlaskTask(Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery_app = Celery(app.import_name)
    celery_app.conf.update(
        broker_url=app.config["CELERY_BROKER_URL"],
        result_backend=app.config["CELERY_RESULT_BACKEND"],
        task_serializer="json",
        result_serializer="json",
        accept_content=["json"],
        task_ignore_result=False,
        worker_hijack_root_logger=False,
        worker_redirect_stdouts=False,
    )
    celery_app.Task = FlaskTask
    celery_app.autodiscover_tasks(["flask_backend.celery_tasks"])
    app.extensions["celery"] = celery_app
    return celery_app
