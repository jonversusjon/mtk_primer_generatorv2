import os
from tests.test_data import TEST_SEQ, TEST_TEMPLATE_SEQ

class Config:
    SECRET_KEY = os.getenv("SECRET_KEY", "default_secret_key")
    DEBUG = True
    TESTING = False

class TestConfig(Config):
    TESTING = True
    TEST_SEQ = TEST_SEQ
    TEST_TEMPLATE_SEQ = TEST_TEMPLATE_SEQ