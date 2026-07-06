"""Entry point: python -m webapp -> serves the Dash app on localhost:8050."""
from webapp.server import build_app

def main() -> None:
    build_app().run(host="127.0.0.1", port=8050, debug=False)

if __name__ == "__main__":
    main()
