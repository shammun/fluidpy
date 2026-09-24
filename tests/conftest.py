"""Make the repository root importable so tests can `from fluidpy import ...` and `from tools import ...`."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def pytest_configure(config):
    """Register the ``slow`` marker (subprocess script runs); deselect with ``-m "not slow"``."""
    config.addinivalue_line("markers", "slow: long-running check (e.g. runs every chapter script); -m 'not slow' skips it")
