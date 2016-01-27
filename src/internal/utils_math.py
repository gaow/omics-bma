#! /usr/bin/env python3
# utils_math.py
# Gao Wang (c) 2015
import numpy as np

class PosteriorCalculator(object):
    def __init__(self, beta):
        pass

def test_almost_equal_recursive(x, y, level = 6, level_cutoff = 3):
  try:
      np.testing.assert_almost_equal(x, y, level)
  except AssertionError:
      if level <= level_cutoff:
          raise
      else:
          level = level - 1
          test_almost_equal_recursive(x, y, level)
  return level
