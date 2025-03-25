import os
import logging

# --- Application and Flask Configurations ---
class BaseConfig:
    SECRET_KEY = os.getenv("SECRET_KEY", "default_secret_key")
    DEBUG = False
    TESTING = False

    # Redis + Celery configurations
    CELERY_BROKER_URL = os.getenv("CELERY_BROKER_URL")
    CELERY_RESULT_BACKEND = os.getenv("CELERY_RESULT_BACKEND")
    CELERY_TASK_SERIALIZER = "json"
    CELERY_RESULT_SERIALIZER = "json"
    CELERY_ACCEPT_CONTENT = ["json"]
    CELERY_IGNORE_RESULT = False
    CELERY_WORKER_HIJACK_ROOT_LOGGER = False
    CELERY_WORKER_REDIRECT_STDOUTS = False
    CELERY_LOG_LEVEL = os.getenv("CELERY_LOG_LEVEL", "INFO")
    CELERY_WORKER_LOG_FILE = os.getenv("CELERY_WORKER_LOG_FILE", "")
    CELERY_LOG_FORMAT = "[%(asctime)s: %(levelname)s/%(processName)s] %(message)s"

    # Redis URL for flask-sse
    REDIS_URL = os.getenv("REDIS_URL", "redis://redis:6379/0")

    # CORS settings
    CORS_ORIGINS = os.getenv("CORS_ORIGINS", "").split(",")
    CORS_METHODS = ["GET", "POST", "PUT", "DELETE", "OPTIONS"]

class DevelopmentConfig(BaseConfig):
    DEBUG = True
    CORS_ORIGINS = ["http://localhost:3000"]

class ProductionConfig(BaseConfig):
    CELERY_WORKER_LOG_FILE = "/var/log/myapp/celery.log"

class TestConfig(BaseConfig):
    TESTING = True
    DEBUG = True

# --- Logging Configuration ---
DEBUG_MODE = os.getenv("DEBUG_MODE", "False").lower() in ("true", "1", "yes")
logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s - %(module)s - %(message)s",
    level=logging.DEBUG if DEBUG_MODE else logging.INFO
)

logger = logging.getLogger("GoldenGateApp")