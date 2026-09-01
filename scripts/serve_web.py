#!/usr/bin/env python3
"""Serve an EDMD web build with WebAssembly-friendly response headers."""

from __future__ import annotations

import argparse
import functools
import http.server
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DIRECTORY = PROJECT_ROOT / "docs"


class WebBuildRequestHandler(http.server.SimpleHTTPRequestHandler):
    """Static-file handler with deterministic web-build content types."""

    def guess_type(self, path: str) -> str:
        lower_path = path.lower()
        if lower_path.endswith(".wasm"):
            return "application/wasm"
        if lower_path.endswith(".js") or lower_path.endswith(".mjs"):
            return "text/javascript"
        if lower_path.endswith(".data"):
            return "application/octet-stream"
        return super().guess_type(path)

    def end_headers(self) -> None:
        # These headers allow a pthread-enabled WebAssembly build to use
        # SharedArrayBuffer while remaining harmless for single-threaded builds.
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        self.send_header("Cross-Origin-Resource-Policy", "same-origin")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Cache-Control", "no-store")
        super().end_headers()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Serve the EDMD browser build (docs/ by default)."
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default=DEFAULT_DIRECTORY,
        type=Path,
        help="directory to serve (default: %(default)s)",
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="interface to bind (default: %(default)s)",
    )
    parser.add_argument(
        "--port",
        default=8000,
        type=int,
        help="TCP port; use 0 to select a free port (default: %(default)s)",
    )
    args = parser.parse_args()
    if not 0 <= args.port <= 65535:
        parser.error("--port must be between 0 and 65535")
    return args


def display_host(host: str) -> str:
    if host in {"", "0.0.0.0", "::"}:
        return "localhost"
    if ":" in host and not host.startswith("["):
        return f"[{host}]"
    return host


def main() -> int:
    args = parse_args()
    directory = args.directory.expanduser().resolve()
    if not directory.is_dir():
        print(f"serve_web.py: not a directory: {directory}", file=sys.stderr)
        return 2

    handler = functools.partial(WebBuildRequestHandler, directory=str(directory))
    try:
        server = http.server.ThreadingHTTPServer((args.host, args.port), handler)
    except OSError as error:
        print(f"serve_web.py: could not start server: {error}", file=sys.stderr)
        return 1

    host, port = server.server_address[:2]
    visible_host = display_host(args.host or str(host))
    print(f"Serving {directory} at http://{visible_host}:{port}/", flush=True)
    print("Press Ctrl-C to stop.", flush=True)

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping server.", flush=True)
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
