#!/usr/bin/env python3
"""Phasor analysis of a LINCam FLIM recording with PhasorPy.

The ``photonscore`` package reads the ``.photons`` file; everything after the
histogram is PhasorPy (https://www.phasorpy.org). For every photon the file
holds

    /photons/x, /photons/y   uint16  position, 12 bit
    /photons/dt              uint16  arrival after the laser pulse, 12 bit,
                                     /photons/TacChannel ps a channel

and the analysis is

    1. histogram the photons into a (y, x, time) cube   photonscore, numpy
    2. the phasor of every pixel's decay                phasor_from_signal
    3. calibrate against the instrument response        phasor_calibrate
    4. drop dim pixels, median-filter the rest          phasor_threshold,
                                                        phasor_filter_median
    5. the populations on the phasor plot               phasor_cluster_gmm
    6. two lifetimes and their share in every pixel     phasor_semicircle_
                                                        intersect, phasor_
                                                        component_fraction

The 4096 TAC channels are taken as one period, 50.4 ns for 12.3 ps channels:
the decay has died away well inside it, so the histogram is one period of a
periodic signal, and its first harmonic is 1 / 50.4 ns = 19.8 MHz.

Usage
-----
    python phasor.py convallaria.photons
    python phasor.py convallaria.photons --harmonic 2 --clusters 4
    python phasor.py convallaria.photons --pixels 1024 --min-photons 100

Writes ``<input>.phasor.png`` with the four panels, and
``<input>.phasor.ome.tif`` with the calibrated phasors of the first three
harmonics, which ``phasorpy.io.phasor_from_ometiff`` reads back. The file is
read in chunks, so the memory in use does not grow with it.
"""

import argparse
import os
import sys
import time

import numpy as np

from photonscore import Vault

HARMONICS = (1, 2, 3)     # computed, calibrated and saved; one is plotted
CHUNK = 16_000_000        # photons per read: ~1.5 GB of temporaries at most


def chunks(count, size=CHUNK):
  """(start, stop) index pairs covering ``range(count)``."""
  for lo in range(0, count, size):
    yield lo, min(lo + size, count)


def geometry(vault):
  """``(positions, channels, channel_ps)`` of the file: 4096, 4096, 12.3."""
  attrs = vault.attributes
  return (1 << int(attrs.get("/photons/PositionBits", 12)),
          1 << int(attrs.get("/photons/TacBits", 12)),
          float(attrs["/photons/TacChannel"]))


def read(vault, pixels=512, bins=256, progress=None):
  """Histogram every photon by pixel and arrival time.

  Returns ``(cube, decay)``: ``cube[y, x, t]``, ``pixels`` a side and
  ``bins`` in time, and the decay of all photons at the full resolution of
  the file, one bin a TAC channel, which is where time zero is read.
  """
  positions, channels, _ = geometry(vault)
  photons = vault.data.photons
  count = photons.x.shape[0]
  cube = np.zeros(pixels * pixels * bins, dtype=np.uint32)
  decay = np.zeros(channels, dtype=np.int64)
  for lo, hi in chunks(count):
    x = photons.x[lo:hi].astype(np.int64)
    y = photons.y[lo:hi].astype(np.int64)
    dt = photons.dt[lo:hi].astype(np.int64)
    voxel = ((y * pixels // positions) * pixels
             + x * pixels // positions) * bins + dt * bins // channels
    np.add(cube, np.bincount(voxel, minlength=cube.size), out=cube,
           casting="unsafe")
    decay += np.bincount(dt, minlength=channels)
    if progress:
      progress(hi, count)
  return cube.reshape(pixels, pixels, bins), decay


def time_zero(decay, channel_ps):
  """Time zero in ns: the steepest rise of the decay, between channels.

  LINCam's response is some 80 ps wide, and its centre is where the decay
  of all photons climbs fastest; on simulated photons this lands within a
  few ps of it for lifetimes from 0.5 to 4 ns.
  """
  from scipy.ndimage import gaussian_filter1d

  slope = np.gradient(gaussian_filter1d(decay.astype(float), 1.0))
  k = int(np.argmax(slope))
  a, b, c = slope[k - 1], slope[k], slope[k + 1]
  return (k + 0.5 * (a - c) / (a - 2 * b + c) + 0.5) * channel_ps * 1e-3


def phasors(cube, period_ns, t0_ns, harmonics=HARMONICS):
  """Calibrated phasors of every pixel: ``(mean, real, imag)``.

  ``mean`` is photons per time bin; ``real`` and ``imag`` have one plane per
  harmonic. The reference is the instrument response: far narrower than a
  period, its phasor is a point on the unit circle at the phase of time
  zero, a reference of lifetime 0.
  """
  from phasorpy.lifetime import phasor_calibrate
  from phasorpy.phasor import phasor_from_polar, phasor_from_signal

  harmonics = list(harmonics)
  frequency = 1e3 / period_ns
  mean, real, imag = phasor_from_signal(cube, axis=-1, harmonic=harmonics)
  # phasor_from_signal puts bin k at k * period / bins, half a bin before
  # the photons in it arrive on average.
  t0 = t0_ns - 0.5 * period_ns / cube.shape[-1]
  phase = 2 * np.pi * np.array(harmonics) * frequency * 1e-3 * t0
  reference = phasor_from_polar(phase, np.ones_like(phase))
  real, imag = phasor_calibrate(real, imag, 1.0, *reference,
                                frequency=frequency, lifetime=0,
                                harmonic=harmonics)
  return mean, real, imag


def clean(mean, real, imag, bins, min_photons=300, repeat=2):
  """NaN where a pixel has fewer photons, then a 3 x 3 median, twice."""
  from phasorpy.filter import phasor_filter_median, phasor_threshold

  mean, real, imag = phasor_threshold(mean, real, imag,
                                      mean_min=min_photons / bins)
  return phasor_filter_median(mean, real, imag, size=3, repeat=repeat)


def clusters(real, imag, count=3, sigma=1.5):
  """Gaussian mixture on the phasor plot, shortest lifetime first.

  Returns the ellipses, ``(centre_real, centre_imag, radius, radius_minor,
  angle)`` at ``sigma``, and a mask per ellipse. A pixel inside two ellipses
  goes to the one it is closer to in units of that ellipse's size, so every
  pixel has at most one colour.
  """
  from phasorpy.cluster import phasor_cluster_gmm
  from phasorpy.cursor import mask_from_elliptic_cursor

  ellipses = phasor_cluster_gmm(real, imag, clusters=count, sigma=sigma,
                                sort="polar", random_state=0)
  masks = np.array([
    mask_from_elliptic_cursor(real, imag, g, s, radius=r, radius_minor=rm,
                              angle=a)
    for g, s, r, rm, a in zip(*ellipses)])
  distance = np.array([
    _elliptic_distance(real, imag, g, s, r, rm, a)
    for g, s, r, rm, a in zip(*ellipses)])
  distance[~masks] = np.inf
  nearest = np.argmin(distance, axis=0)
  masks &= nearest == np.arange(count)[:, None, None]
  return ellipses, masks


def _elliptic_distance(real, imag, g, s, radius, radius_minor, angle):
  ca, sa = np.cos(angle), np.sin(angle)
  u, v = real - g, imag - s
  return np.hypot((u * ca + v * sa) / radius, (v * ca - u * sa) / radius_minor)


def components(mean, real, imag, frequency):
  """Two lifetimes whose mixtures lie along the phasor cloud.

  The line is the cloud's long axis, photon-weighted; where it meets the
  semicircle are the phasors of two single exponentials. Returns
  ``(lifetimes, fraction)``: the two lifetimes in ns, shorter first, and per
  pixel the share of the photons from the shorter one.
  """
  from phasorpy.component import phasor_component_fraction
  from phasorpy.lifetime import (phasor_from_lifetime,
                                 phasor_semicircle_intersect,
                                 phasor_to_apparent_lifetime)

  ok = np.isfinite(real)
  w = mean[ok]
  g0 = np.average(real[ok], weights=w)
  s0 = np.average(imag[ok], weights=w)
  cov = np.cov(np.vstack((real[ok] - g0, imag[ok] - s0)), aweights=w)
  dg, ds = np.linalg.eigh(cov)[1][:, -1]
  ends = phasor_semicircle_intersect(g0 - dg, s0 - ds, g0 + dg, s0 + ds)
  lifetimes = sorted(
    float(phasor_to_apparent_lifetime(ends[i], ends[i + 1], frequency)[0])
    for i in (0, 2))
  fraction = phasor_component_fraction(
    real, imag, *phasor_from_lifetime(frequency, lifetimes))
  return lifetimes, fraction


# --- The figure, shared with the notebook -----------------------------------

COLORS = ("#3a7bd5", "#2ca25f", "#CC3333", "#e6a23c", "#8e6cc9")


def preview(image):
  """As Photonscore Preview shows an image: x from right to left."""
  return image[:, ::-1]


def brightness(counts, gamma=0.5, percentile=99.9):
  """Counts as 0 to 1, square-root scaled so dim tissue still shows."""
  white = max(1.0, float(np.percentile(counts, percentile)))
  return np.clip(counts / white, 0, 1) ** gamma


def intensity_panel(ax, counts, title=None):
  ax.imshow(preview(brightness(counts)), cmap="gray", vmin=0, vmax=1,
            interpolation="none")
  ax.set_title(title or f"{int(counts.sum()):,} photons")
  ax.set_axis_off()


def phasor_panel(ax, real, imag, frequency, ellipses=None, lifetimes=None,
                 harmonic=None, dark=False):
  from phasorpy.lifetime import phasor_from_lifetime
  from phasorpy.plot import PhasorPlot

  from matplotlib.colors import LinearSegmentedColormap

  ink = "#d0d0d0" if dark else "black"
  # From a grey the page does not swallow, so one pixel in a bin still shows.
  shades = ("#5a5a5a", "#ffffff") if dark else ("#b4b4b4", "#000000")
  label = f"harmonic {harmonic}, " if harmonic else ""
  plot = PhasorPlot(ax=ax, frequency=frequency,
                    title=f"Phasor plot, {label}{frequency:.1f} MHz")
  plot.hist2d(real, imag, bins=300,
              cmap=LinearSegmentedColormap.from_list("counts", shades))
  plot.semicircle(frequency=frequency, lifetime=[0.25, 0.5, 1, 2, 4, 8])
  if ellipses is not None:
    for i, (g, s, r, rm, a) in enumerate(zip(*ellipses)):
      plot.cursor(g, s, radius=r, radius_minor=rm, angle=a,
                  color=COLORS[i % len(COLORS)], linewidth=1.5)
  if lifetimes is not None:
    g, s = phasor_from_lifetime(frequency, lifetimes)
    ax.plot(g, s, "-o", color=ink, lw=1, ms=4)
    for gi, si, tau in zip(g, s, lifetimes):
      ax.annotate(f"{tau:.2f} ns", (gi, si), color=ink,
                  textcoords="offset points", xytext=(6, 6))
  return plot


def clusters_panel(ax, counts, masks):
  """Each ellipse's pixels in its colour, the rest in grey."""
  from matplotlib.colors import to_rgb
  from phasorpy.cursor import pseudo_color

  lum = brightness(counts)
  colors = [to_rgb(COLORS[i % len(COLORS)]) for i in range(len(masks))]
  rgb = pseudo_color(*masks, intensity=0.3 + 0.7 * lum, colors=colors,
                     vmin=0, vmax=1)
  rest = ~masks.any(axis=0)
  rgb[rest] = 0.6 * lum[rest, None]
  ax.imshow(preview(rgb), interpolation="none")
  ax.set_title("Pixels in each ellipse")
  ax.set_axis_off()


def fraction_panel(ax, counts, fraction, lifetimes, cmap="coolwarm"):
  """The share as colour, the photons as brightness."""
  import matplotlib as mpl

  lum = brightness(counts)
  rgb = mpl.colormaps[cmap](np.nan_to_num(fraction, nan=0.5))[..., :3]
  rgb *= (0.3 + 0.7 * lum)[..., None]
  rest = ~np.isfinite(fraction)
  rgb[rest] = 0.6 * lum[rest, None]
  ax.imshow(preview(rgb), interpolation="none")
  ax.set_title(f"Share of the {lifetimes[0]:.2f} ns component")
  ax.set_axis_off()
  scale = mpl.cm.ScalarMappable(mpl.colors.Normalize(0, 1), cmap)
  ax.figure.colorbar(scale, ax=ax, shrink=0.7)


def figure(counts, real, imag, frequency, ellipses, masks, lifetimes,
           fraction, harmonic=None, dark=False):
  """The four panels of the note: image, phasor plot, clusters, share."""
  import matplotlib.pyplot as plt

  with plt.style.context("dark_background" if dark else "default"):
    fig, axes = plt.subplots(2, 2, figsize=(12, 12), layout="constrained")
    intensity_panel(axes[0, 0], counts)
    phasor_panel(axes[0, 1], real, imag, frequency, ellipses, lifetimes,
                 harmonic, dark)
    clusters_panel(axes[1, 0], counts, masks)
    fraction_panel(axes[1, 1], counts, fraction, lifetimes)
  return fig


def _progress(label):
  began = time.time()

  def show(done, count):
    print(f"\r  {label}: {100 * done / count:5.1f} %"
          f"  {time.time() - began:5.1f} s",
          end="\n" if done == count else "", flush=True)
  return show


def main(argv=None):
  parser = argparse.ArgumentParser(
    description="Phasor analysis of a LINCam .photons file with PhasorPy.")
  parser.add_argument("input", help="the .photons file")
  parser.add_argument("--pixels", type=int, default=512,
                      help="pixels along each side of the image (512)")
  parser.add_argument("--bins", type=int, default=256,
                      help="time bins over the 4096 TAC channels (256)")
  parser.add_argument("--harmonic", type=int, default=3, choices=HARMONICS,
                      help="the harmonic plotted and analysed (3)")
  parser.add_argument("--min-photons", type=int, default=300, metavar="N",
                      help="pixels with fewer photons are left out (300)")
  parser.add_argument("--clusters", type=int, default=3,
                      help="populations sought on the phasor plot (3)")
  parser.add_argument("--dark", action="store_true",
                      help="draw the figure for a dark page")
  args = parser.parse_args(argv)

  vault = Vault(args.input)
  positions, channels, channel_ps = geometry(vault)
  period_ns = channels * channel_ps * 1e-3
  print(f"{args.input}: {vault.data.photons.x.shape[0]:,} photons, "
        f"{channel_ps} ps a channel, a period of {period_ns:.2f} ns")

  cube, decay = read(vault, args.pixels, args.bins, _progress("read"))
  vault.close()
  t0 = time_zero(decay, channel_ps)
  print(f"  time zero at {t0:.3f} ns")

  mean, real, imag = phasors(cube, period_ns, t0)
  counts = cube.sum(axis=-1)
  del cube
  mean, real, imag = clean(mean, real, imag, args.bins, args.min_photons)
  h = HARMONICS.index(args.harmonic)
  frequency = args.harmonic * 1e3 / period_ns
  print(f"  {int(np.isfinite(real[h]).sum()):,} pixels with "
        f"{args.min_photons} photons or more; harmonic {args.harmonic} "
        f"is {frequency:.1f} MHz")

  from phasorpy.lifetime import phasor_to_apparent_lifetime

  ellipses, masks = clusters(real[h], imag[h], args.clusters)
  for i, (g, s) in enumerate(zip(ellipses[0], ellipses[1])):
    phase, modulation = phasor_to_apparent_lifetime(g, s, frequency)
    print(f"  cluster {i + 1}: G {g:.3f}, S {s:.3f}, tau_phase "
          f"{float(phase):.2f} ns, tau_mod {float(modulation):.2f} ns, "
          f"{int(masks[i].sum()):,} pixels")
  lifetimes, fraction = components(mean, real[h], imag[h], frequency)
  print(f"  two components: {lifetimes[0]:.2f} ns and {lifetimes[1]:.2f} ns")

  stem = os.path.splitext(args.input)[0]
  fig = figure(counts, real[h], imag[h], frequency, ellipses, masks,
               lifetimes, fraction, args.harmonic, args.dark)
  png = stem + (".phasor-dark.png" if args.dark else ".phasor.png")
  fig.savefig(png, dpi=150, transparent=True)
  print(f"  figure in {png}")

  from phasorpy.io import phasor_to_ometiff

  tif = stem + ".phasor.ome.tif"
  phasor_to_ometiff(tif, mean, real, imag, frequency=1e3 / period_ns,
                    harmonic=list(HARMONICS), dtype=np.float32,
                    description=f"{os.path.basename(args.input)}, "
                                f"time zero {t0:.4f} ns")
  print(f"  phasors in {tif}")


if __name__ == "__main__":
  sys.exit(main())
