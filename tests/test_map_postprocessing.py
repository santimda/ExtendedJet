"""Fast checks for the map post-processing conventions."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "postprocessing"))
import Fmap  # noqa: E402


class MapPostprocessingTests(unittest.TestCase):
    def setUp(self):
        self.pixel_mas = Fmap.PIXEL_MAS
        self.beam = (0.31, 0.31)
        self.shape = (101, 101)
        self.pixel_intensity = 1e-26 / (self.pixel_mas * Fmap.MAS) ** 2

    def test_single_pixel_is_one_mJy_per_beam(self):
        image = np.zeros(self.shape)
        image[50, 50] = self.pixel_intensity
        convolved = Fmap.convolve(image, self.beam, self.pixel_mas)
        self.assertAlmostEqual(float(convolved.max()), 1.0, places=12)

    def test_nearby_pixels_add_in_mJy_per_beam(self):
        image = np.zeros(self.shape)
        image[50, 50] = self.pixel_intensity
        image[50, 51] = self.pixel_intensity
        convolved = Fmap.convolve(image, self.beam, self.pixel_mas)
        self.assertGreater(float(convolved.max()), 1.99)
        self.assertLess(float(convolved.max()), 2.0)


if __name__ == "__main__":
    unittest.main()
