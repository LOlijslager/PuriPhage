"""
logging_utils.py

Shared logging configuration used throughout the pipeline.
"""

import logging
import os


def setup_logging(log_dir, sample_name, troubleshooting=False):
    """
    Configure console and file logging.

    Returns
    -------
    logger
    """

    os.makedirs(log_dir, exist_ok=True)

    logfile = os.path.join(
        log_dir,
        f"{sample_name}.log"
    )

    if os.path.exists(logfile):
        os.remove(logfile)

    level = (
        logging.DEBUG
        if troubleshooting
        else logging.INFO
    )

    logging.basicConfig(
        level=level,
        format=(
            "%(asctime)s - "
            "%(levelname)s - "
            "%(message)s"
        ),
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(logfile),
            logging.StreamHandler(),
        ],
        force=True,
    )

    return logging.getLogger("phage_pipeline")
