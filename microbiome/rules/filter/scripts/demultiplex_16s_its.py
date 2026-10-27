#!/usr/bin/env python
"""Compatibility entry point for the shared QIAseq region demultiplexer."""

from pathlib import Path
import runpy


if __name__ == '__main__':
    runpy.run_path(str(Path(__file__).resolve().parents[4] / 'common/scripts/demultiplex_16s_its.py'), run_name='__main__')
