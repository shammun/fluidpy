"""Make the repository root importable so tests can `from fluidpy import ...` and `from tools import ...`."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
