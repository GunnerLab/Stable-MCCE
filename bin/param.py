#!/usr/bin/env python
"""
This module handles MCCE topology parameters and run-time control parameters
"""

import logging
logging.basicConfig(level=logging.INFO, format='%(levelname)-s: %(message)s')

class Env:
    def __init__(self):
        # hard code values
        self.tpl = {}
        return

    def load_ftpl(self, fname):
        """Load a tpl file"""



