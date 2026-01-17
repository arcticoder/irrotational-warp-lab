from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np


def _to_jsonable(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return {
            "shape": list(obj.shape),
            "min": float(np.min(obj)),
            "max": float(np.max(obj)),
            "mean": float(np.mean(obj)),
        }
    return obj


def write_summary_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cooked = {k: _to_jsonable(v) for k, v in payload.items()}
    path.write_text(json.dumps(cooked, indent=2, sort_keys=True) + "\n")
