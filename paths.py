from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent

SRC_DIR = PROJECT_ROOT / "src"
RESULTS_DIR = PROJECT_ROOT / "results"

DIST_DIR = RESULTS_DIR / "distribution"
PARA_DIR = RESULTS_DIR / "parameters"
ATK_DIR = RESULTS_DIR / "attack"
DATA_DIR = RESULTS_DIR / "data"

DIST_DIR.mkdir(parents=True, exist_ok=True)
PARA_DIR.mkdir(parents=True, exist_ok=True)
ATK_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)
