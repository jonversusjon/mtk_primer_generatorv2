import logging
import os

# Get debug mode from environment variable
DEBUG_MODE = os.getenv("DEBUG_MODE", "False").lower() in ("true", "1", "yes")

# Configure a single logging instance
logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
    level=logging.DEBUG if DEBUG_MODE else logging.INFO
)

# Create a centralized logger
logger = logging.getLogger("GoldenGateApp")
