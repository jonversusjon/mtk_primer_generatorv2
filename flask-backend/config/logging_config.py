import logging
import os

# Get debug mode from environment variable
DEBUG_MODE = os.getenv("DEBUG_MODE", "False").lower() in ("true", "1", "yes")

# Configure a single logging instance with a format that includes module info.
logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s - %(module)s - %(message)s",
    level=logging.DEBUG if DEBUG_MODE else logging.INFO
)

# Create a centralized logger for the app.
logger = logging.getLogger("GoldenGateApp")
