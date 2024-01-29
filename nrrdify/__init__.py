#!/usr/bin/env python

# ========================================================================
#  Copyright Het Nederlands Kanker Instituut - Antoni van Leeuwenhoek
#
#  Licensed under the 3-clause BSD License
# ========================================================================

import logging

logger = logging.getLogger('nrrdify')
_post_processing = None


def setPostProcessing(func):
    global _post_processing
    _post_processing = func


def getPostProcessing():
    global _post_processing
    return _post_processing


if len(logger.handlers) == 0 and len(logging.getLogger().handlers) == 0:
  print('Adding handler for logger')
  handler = logging.StreamHandler()
  formatter = logging.Formatter('[%(asctime)-.19s] %(levelname)-.1s: %(message)s')
  handler.setFormatter(formatter)
  handler.setLevel(logging.INFO)

  logger.addHandler(handler)
  logger.setLevel(logging.INFO)
