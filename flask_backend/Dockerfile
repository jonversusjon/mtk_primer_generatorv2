# Use a lightweight Python base
FROM python:3.10-slim

# Prevent Python buffering (helps logs)
ENV PYTHONUNBUFFERED=1

# Set the working directory to the repository root
WORKDIR /app

# Copy the requirements file from flask_backend for dependency caching
COPY flask_backend/requirements.txt flask_backend/requirements.txt
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r flask_backend/requirements.txt
    
# Install dependencies
RUN pip install --no-cache-dir -r flask_backend/requirements.txt

# Copy the entire repository into the container
COPY . /app

# Expose the Flask port
EXPOSE 5000

# Default command (Flask dev server)
CMD ["flask", "run", "--host=0.0.0.0", "--reload"]