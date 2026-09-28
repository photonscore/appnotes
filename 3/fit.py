#!/usr/bin/env python3
"""Lifetime fitting of a LINCam FLIM recording with FLIMKit.

FLIMKit (https://github.com/FLIMKit/FLIMKit) opens the ``.photons`` file
itself, through ``photonsfile``, a pure-Python reader of the D7 format, and
bins every photon by position and arrival time. From there it is FLIMKit's
reconvolution fit:

    1. read the photons into a (y, x, time) cube      FLIMFile
    2. the instrument response, a Gaussian            gaussian_irf_from_fwhm
    3. fit the decay of all photons from its rise     fit_summed
    4. fit every pixel, the lifetimes held            fit_per_pixel
    5. the maps: mean lifetime, and the share of each lifetime

The fit starts where the decay of all photons climbs through half its peak.
The rise before it is where LINCam's response is no Gaussian: a core some
70 ps wide behind a weak shoulder, a thousandth of the peak and 1 ns long.

Usage
-----
    python fit.py convallaria.photons
    python fit.py convallaria.photons --nexp 3 --fit-start 10.5
    python fit.py convallaria.photons --binning 1 --min-photons 250

Writes ``<input>.fit.png`` with the four panels, and ``<input>.fit.tif`` with
the maps as 32-bit float planes. FLIMKit reads the whole file into memory: a
recording of 284 million photons takes some 24 GB.
"""

import argparse
import os
import sys
import time

import numpy as np


def load(path, pixels=512, bins=1024):
  """The file as FLIMKit opens it: ``(flim, decay)``.

  ``flim`` is FLIMKit's reader, which ``pixel_stack`` and ``intensity_image``
  ask for the cube and the image; ``decay`` is the decay of all photons.
  """
  from flimkit.formats import FLIMFile

  flim = FLIMFile(path, verbose=True, pixels=pixels, n_bins=bins)
  return flim, flim.summed_decay()


def response(decay, bin_s, fwhm_ns=0.08):
  """A Gaussian instrument response, ``fwhm_ns`` wide, at the decay's peak.

  FLIMKit shifts and widens it in the fit of the decay of all photons; the
  80 ps is LINCam's own, on the scale of a 49 ps bin.
  """
  from flimkit.FLIM.irf_tools import gaussian_irf_from_fwhm

  return gaussian_irf_from_fwhm(decay.size, bin_s, fwhm_ns,
                                int(np.argmax(decay)))


def fit_start(decay, bin_s, level=0.5):
  """Where the decay first reaches ``level`` of its peak, in ns."""
  return int(np.argmax(decay >= level * decay.max())) * bin_s * 1e9


def fit_all(decay, bin_s, irf, nexp=4, start_ns=None, tau_min_ns=0.05,
            tau_max_ns=20.0):
  """Reconvolution fit of the decay of all photons: ``(popt, summary)``.

  Differential evolution, then Levenberg-Marquardt, on the Poisson
  likelihood; the response's shift and width and a flat background are
  fitted with the lifetimes, from ``start_ns`` on.
  """
  from flimkit.FLIM.fitters import fit_summed

  return fit_summed(decay, bin_s, decay.size, irf, False, True, True, nexp,
                    tau_min_ns, tau_max_ns, optimizer="de", workers=1,
                    fit_start_ns=start_ns)


def fit_pixels(stack, bin_s, irf, popt, summary, nexp=4, min_photons=1000):
  """Every pixel with the lifetimes of the decay of all photons held.

  Only the amplitudes are free, one per lifetime: FLIMKit solves each pixel
  by non-negative least squares over the same fit window.
  """
  from flimkit.FLIM.fitters import fit_per_pixel

  return fit_per_pixel(stack, bin_s, stack.shape[-1], irf, False, True, True,
                       popt, nexp, min_photons=min_photons,
                       fit_idx=summary["fit_idx"])


def shares(maps, taus_ns):
  """Per pixel, the share of the photons from each lifetime.

  FLIMKit's ``frac_i`` are amplitude fractions; the photons of component i
  are its amplitude times its lifetime.
  """
  photons = np.array([maps[f"alpha_{i + 1}"] * tau
                      for i, tau in enumerate(taus_ns)])
  with np.errstate(invalid="ignore", divide="ignore"):
    return photons / photons.sum(axis=0)


# --- The figure, shared with the notebook -----------------------------------

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


def decay_panel(ax, below, decay, bin_s, summary, dark=False):
  """The decay of all photons and the fit over its window; the residuals,
  in standard deviations, on ``below``."""
  ink = "#d0d0d0" if dark else "black"
  t = (np.arange(decay.size) + 0.5) * bin_s * 1e9
  window = np.asarray(summary["fit_idx"])
  model = np.asarray(summary["model"])[window]
  taus = ", ".join(f"{tau:.2f}" for tau in sorted(summary["taus_ns"]))
  ax.semilogy(t, np.maximum(decay, 0.5), ".", ms=2, color="#7f7f7f",
              label="photons")
  ax.semilogy(t[window], model, color="#CC3333", lw=1.2,
              label=f"fit: {taus} ns")
  ax.set_xlim(t[window[0]] - 2, t[window[-1]] + 1)
  ax.set_ylim(max(1.0, decay[window].min() / 3), decay.max() * 2)
  ax.set_ylabel(f"photons per {bin_s * 1e12:.0f} ps")
  ax.set_title("All photons")
  ax.legend(loc="upper right", frameon=False, labelcolor=ink)
  ax.tick_params(labelbottom=False)
  residual = (decay[window] - model) / np.sqrt(np.maximum(model, 1))
  below.plot(t[window], residual, lw=0.6, color=ink)
  below.axhline(0, color="#CC3333", lw=0.8)
  # The rise, where the response is no Gaussian, would set the scale alone.
  span = max(5.0, float(np.percentile(np.abs(residual), 99)))
  below.set_ylim(-span, span)
  below.set_ylabel("residual, σ")
  below.set_xlabel("arrival after the laser pulse, ns")


def map_panel(ax, counts, values, title, vmin, vmax, cmap, label):
  """A map as colour, the photons as brightness."""
  import matplotlib as mpl

  lum = brightness(counts)
  scaled = np.clip((np.nan_to_num(values, nan=vmin) - vmin) / (vmax - vmin),
                   0, 1)
  rgb = mpl.colormaps[cmap](scaled)[..., :3] * (0.3 + 0.7 * lum)[..., None]
  rest = ~np.isfinite(values)
  rgb[rest] = 0.6 * lum[rest, None]
  ax.imshow(preview(rgb), interpolation="none")
  ax.set_title(title)
  ax.set_axis_off()
  scale = mpl.cm.ScalarMappable(mpl.colors.Normalize(vmin, vmax), cmap)
  ax.figure.colorbar(scale, ax=ax, shrink=0.7, label=label)


def histogram_panel(ax, counts, values, vmin, vmax, cmap, bins=120):
  """How the photons spread over the map's values, coloured as the map."""
  import matplotlib as mpl

  ok = np.isfinite(values)
  weights, edges = np.histogram(values[ok], bins=bins, range=(vmin, vmax),
                                weights=counts[ok])
  centres = 0.5 * (edges[:-1] + edges[1:])
  colors = mpl.colormaps[cmap]((centres - vmin) / (vmax - vmin))
  ax.bar(centres, weights / weights.sum(), width=edges[1] - edges[0],
         color=colors)
  ax.set_xlim(vmin, vmax)
  ax.set_xlabel("mean lifetime, intensity-weighted, ns")
  ax.set_ylabel("share of the photons")
  ax.set_title("Photons by mean lifetime")


def figure(counts, decay, bin_s, summary, maps, dark=False):
  """The panels of the note: image, decay, mean lifetime, histogram."""
  import matplotlib.pyplot as plt

  tau = maps["tau_mean_int"]
  lo, hi = np.nanpercentile(tau, [1, 99])
  with plt.style.context("dark_background" if dark else "default"):
    fig, axes = plt.subplot_mosaic(
      [["image", "decay"], ["image", "residual"], ["map", "histogram"]],
      height_ratios=[3, 1, 4], figsize=(12, 12), layout="constrained")
    axes["residual"].sharex(axes["decay"])
    intensity_panel(axes["image"], counts)
    decay_panel(axes["decay"], axes["residual"], decay, bin_s, summary, dark)
    map_panel(axes["map"], maps["intensity"], tau,
              "Mean lifetime, intensity-weighted", lo, hi, "turbo", "ns")
    histogram_panel(axes["histogram"], maps["intensity"], tau, lo, hi,
                    "turbo")
  return fig


def save_maps(path, maps, share, taus_ns):
  """The maps as float32 planes of one TIFF, named in its description."""
  import tifffile

  names = ["intensity", "tau_mean_int", "tau_mean_amp", "chi2_r"]
  planes = [maps[n] for n in names]
  for i, tau in enumerate(taus_ns):
    names += [f"alpha_{i + 1}", f"share_{i + 1} ({tau:.3f} ns)"]
    planes += [maps[f"alpha_{i + 1}"], share[i]]
  tifffile.imwrite(path, np.array(planes, dtype=np.float32),
                   photometric="minisblack",
                   description="; ".join(names))
  return names


def main(argv=None):
  parser = argparse.ArgumentParser(
    description="Lifetime fitting of a LINCam .photons file with FLIMKit.")
  parser.add_argument("input", help="the .photons file")
  parser.add_argument("--pixels", type=int, default=512,
                      help="pixels along each side of the image (512)")
  parser.add_argument("--binning", type=int, default=2,
                      help="pixels summed along each side for the fit (2)")
  parser.add_argument("--bins", type=int, default=1024,
                      help="time bins over the 4096 TAC channels (1024)")
  parser.add_argument("--nexp", type=int, default=4, choices=(1, 2, 3, 4),
                      help="exponentials in the model (4)")
  parser.add_argument("--irf-fwhm", type=float, default=0.08, metavar="NS",
                      help="width of the Gaussian response, ns (0.08)")
  parser.add_argument("--fit-start", type=float, metavar="NS",
                      help="where the fit starts, ns (half the peak)")
  parser.add_argument("--min-photons", type=int, default=1000, metavar="N",
                      help="pixels with fewer photons are not fitted (1000)")
  parser.add_argument("--dark", action="store_true",
                      help="draw the figure for a dark page")
  args = parser.parse_args(argv)

  began = time.time()
  flim, decay = load(args.input, args.pixels, args.bins)
  bin_s = flim.tcspc_res
  print(f"  read in {time.time() - began:.0f} s: {int(decay.sum()):,} "
        f"photons, {bin_s * 1e12:.1f} ps a bin")

  irf = response(decay, bin_s, args.irf_fwhm)
  start = args.fit_start or fit_start(decay, bin_s)
  popt, summary = fit_all(decay, bin_s, irf, args.nexp, start)
  taus = np.asarray(summary["taus_ns"])
  print(f"  fit from {start:.2f} ns: lifetimes {np.round(taus, 3)} ns, "
        f"amplitude "
        f"fractions {np.round(summary['fractions'], 3)}, mean lifetime "
        f"{summary['tau_mean_int_ns']:.3f} ns (intensity-weighted), "
        f"reduced chi2 {summary['reduced_chi2']:.1f}")

  counts = flim.intensity_image().astype(float)
  stack = flim.pixel_stack(binning=args.binning)
  began = time.time()
  maps = fit_pixels(stack, bin_s, irf, popt, summary, args.nexp,
                    args.min_photons)
  fitted = np.isfinite(maps["tau_mean_int"])
  share = shares(maps, taus)
  print(f"  {int(fitted.sum()):,} pixels of {stack.shape[0]} x "
        f"{stack.shape[1]} fitted in {time.time() - began:.0f} s, median "
        f"reduced chi2 {np.nanmedian(maps['chi2_r']):.2f}")

  stem = os.path.splitext(args.input)[0]
  fig = figure(counts, decay, bin_s, summary, maps, args.dark)
  png = stem + (".fit-dark.png" if args.dark else ".fit.png")
  fig.savefig(png, dpi=150, transparent=True)
  print(f"  figure in {png}")
  tif = stem + ".fit.tif"
  names = save_maps(tif, maps, share, taus)
  print(f"  maps in {tif}: {', '.join(names)}")


if __name__ == "__main__":
  sys.exit(main())
