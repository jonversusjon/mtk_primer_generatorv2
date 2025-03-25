# custom_sse.py
from flask import Blueprint, Response, request
import redis
import os

custom_sse = Blueprint("custom_sse", __name__)

# Connect to Redis using the same config as your Flask app / Celery worker
r = redis.Redis.from_url(os.environ.get("REDIS_URL", "redis://localhost:6379"))

@custom_sse.route("/stream")
def stream():
    channel = request.args.get("channel", "default")
    pubsub = r.pubsub()
    pubsub.subscribe(channel)

    def event_stream():
        for message in pubsub.listen():
            if message["type"] == "message":
                yield f"data: {message['data'].decode()}\n\n"

    return Response(event_stream(), mimetype="text/event-stream")
