from flask import Flask
from routes.main import main
from routes.api import api_blueprint

app = Flask(__name__)
app.config.from_object("config.settings")

# Register Blueprints
app.register_blueprint(main)
app.register_blueprint(api_blueprint, url_prefix="/api")

if __name__ == "__main__":
    app.run(debug=True)
