# -*- coding: utf-8 -*-
"""Resolve isolated LineaLinkR 0.6.4 work and read-only input roots."""
from __future__ import annotations

import os
from pathlib import Path

CONFIG = {
    "paths": {
        "work": r"F:\biotrade\0.64LLR",
        "input_data": r"F:\biotrade\LineaLinkR_benchmark_assets\01_datasets",
        "baseline_results": r"F:\biotrade\0.64LLR\results",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
}


def work_dir() -> Path:
    value = os.environ.get("LLR064_WORK_ROOT", CONFIG["paths"]["work"])
    return Path(value).expanduser().resolve()


def assets_dir() -> Path:
    return work_dir()


def pipeline_dir() -> Path:
    env = os.environ.get("LLR_PIPELINE", "").strip()
    if env:
        return Path(env).expanduser().resolve()
    return work_dir() / "pipeline"


def input_data_dir() -> Path:
    value = os.environ.get("LLR064_INPUT_DATA", CONFIG["paths"]["input_data"])
    return Path(value).expanduser().resolve()


def baseline_results_dir() -> Path:
    return Path(CONFIG["paths"]["baseline_results"]).resolve()
