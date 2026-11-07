# -*- coding: utf-8 -*-
"""Resolve the current revision work root and read-only benchmark inputs."""
from __future__ import annotations

import os
from pathlib import Path

CONFIG = {
    "paths": {
        "work": str(Path(__file__).resolve().parents[2]),
        "input_data": str(Path(__file__).resolve().parents[2] / "inputs" / "benchmark"),
        "baseline_results": str(Path(__file__).resolve().parents[2] / "inputs" / "comparison" / "results"),
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
}


def work_dir() -> Path:
    value = os.environ.get("LLR_WORK_ROOT", CONFIG["paths"]["work"])
    return Path(value).expanduser().resolve()


def assets_dir() -> Path:
    return work_dir()


def pipeline_dir() -> Path:
    env = os.environ.get("LLR_PIPELINE", "").strip()
    if env:
        return Path(env).expanduser().resolve()
    return work_dir() / "experiments"


def input_data_dir() -> Path:
    value = os.environ.get("LLR_INPUT_DATA", CONFIG["paths"]["input_data"])
    return Path(value).expanduser().resolve()


def baseline_results_dir() -> Path:
    return Path(CONFIG["paths"]["baseline_results"]).resolve()
