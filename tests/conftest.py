import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
for path in (ROOT, ROOT / "bin"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
