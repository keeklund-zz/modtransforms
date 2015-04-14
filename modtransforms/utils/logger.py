"""Module contains logger utility functions."""
import logging
import sys

def build_logger():
    """

    """
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(funcName)s - %(message)s",
        level=logging.INFO)
    return logging
