#!/usr/bin/env python3
from __future__ import annotations
import importlib.util
from pathlib import Path
import common
base = common.RUN_ROOT / "scripts" / "03_integrate_v7_cnv.py"
spec = importlib.util.spec_from_file_location("v7_integrator_relaxed", base)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
module.MANIFEST_ROOT = common.RUN_ROOT / "attempts" / "manifest_v7_relaxed_candidate" / "manifests"
module.main()

