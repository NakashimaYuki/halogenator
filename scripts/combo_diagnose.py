#!/usr/bin/env python
# -*- coding: ascii -*-
"""
CLI wrapper for combinatorial diagnosis module.

Usage:
    python scripts/combo_diagnose.py SMILES --config CONFIG_FILE
"""

import sys
import os

# Add src directory to path so we can import halogenator modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from halogenator.combo_diagnose import main

if __name__ == "__main__":
    sys.exit(main())