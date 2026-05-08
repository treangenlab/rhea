"""
rhea: reference-free heterogeneity and evolution in assembly graphs

Rhea is a software used to detect structural variants (SVs) between steps 
in long-read metagenomic series data.
"""

__author__ = 'Kristen Curry'
__version__ = '1.0.0'
__date__ = 'Jan 2024'

import os
import importlib.util

def main():
    """Entry point for the rhea command-line interface."""
    # Load rhea.py from the parent directory
    rhea_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    rhea_py_path = os.path.join(rhea_dir, 'rhea.py')
    
    spec = importlib.util.spec_from_file_location("rhea_module", rhea_py_path)
    rhea_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rhea_module)
    
    return rhea_module.main()
