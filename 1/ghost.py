#!/usr/bin/env python3
"""Extract a quantum ghost image from a LINCam ``.photons`` file.

The setup is the one of Ryan et al., Optica 11, 1261 (2024): photon pairs from
a PPKTP crystal, the visible signal photon imaged onto the LINCam, the infrared
idler sent through the sample onto a single-pixel "bucket" SPAD. The SPAD's
pulse drives the LINCam's reference input, so the file holds, for every photon
the camera saw:

    /photons/x, /photons/y   uint16  position, 12 bit
    /start/time              uint64  the photon's own arrival time, ps
    /stop/time               uint64  the last bucket click before it, ps
    /photons/ms              uint64  index of the first photon of every ms

The recorded ``/photons/dt`` is the same delay in 100 ps channels, too coarse
for a peak some 150 ps wide, and is not read: every delay here is the TDC's,
``start - stop``, to the picosecond.

The camera image alone shows the bare crystal face: the signal photons never
met the sample. A signal photon whose idler twin made it through the sample
arrives a fixed delay after the bucket click, so ``start - stop`` has a sharp
coincidence peak over a flat floor of accidentals, and the photons under that
peak draw the sample: the ghost image.

    1. histogram every photon                   -> the full image
    2. histogram start - stop                   -> the coincidence peak
    3. keep the photons within +/-100 ps of it  -> the ghost image

Usage
-----
    python ghost.py ghost-dataset.photons
    python ghost.py ghost-dataset.photons --window 150 --save ghost-only.photons
    python ghost.py ghost-dataset.photons --bins 2048 --png

Writes ``<input>.ghost.html`` with the three plots; ``--save`` also writes the
selected photons as a ``.photons`` file that Photonscore Preview, and this
script, open like any other. Everything streams in chunks, so the memory in use
does not grow with the file.
"""

import argparse
import os
import struct
import sys
import time
import zlib
from datetime import datetime

import numpy as np

from photonscore import Vault, histogram

POSITION_MAX = 4096       # x and y are 12 bit
RESOLUTIONS = (512, 1024, 2048)  # pixels a side: 8, 4 or 2 position units each
SPAN_PS = 20_000          # the delay histogram covers [0, SPAN_PS), 1 ps a bin
SMOOTH_PS = 25            # boxcar the peak is searched under
CHUNK = 16_000_000        # photons per read: ~350 MB at its widest


def chunks(count, size=CHUNK):
  """(start, stop) index pairs covering ``range(count)``."""
  for lo in range(0, count, size):
    yield lo, min(lo + size, count)


def delays(vault, lo, hi):
  """``start - stop`` in ps for the photons ``[lo, hi)``, signed."""
  start = vault.data.start.time[lo:hi].astype(np.int64)
  stop = vault.data.stop.time[lo:hi].astype(np.int64)
  return start - stop


def image_of(x, y, bins):
  """Counts as ``[row = y, column = x]``, the way an image is indexed."""
  # Square on purpose: the package fills its 2-D result y-major, which reads as
  # [y, x] only while both axes have the same number of bins.
  return histogram(x, 0, POSITION_MAX, bins, y, 0, POSITION_MAX, bins)


def scan(vault, bins=512, span_ps=SPAN_PS, progress=None):
  """First pass: the image of every photon and the delay histogram.

  Returns ``(full_image, delay_counts)``; ``delay_counts[k]`` is the number of
  photons that came ``k`` ps after the bucket click before them.
  """
  count = vault.data.photons.x.shape[0]
  full = np.zeros((bins, bins), dtype=np.int64)
  delay_counts = np.zeros(span_ps, dtype=np.int64)
  for lo, hi in chunks(count):
    full += image_of(vault.data.photons.x[lo:hi], vault.data.photons.y[lo:hi], bins)
    delay_counts += histogram(delays(vault, lo, hi), 0, span_ps, span_ps)
    if progress:
      progress(hi, count)
  return full, delay_counts


def find_peak(delay_counts, smooth_ps=SMOOTH_PS):
  """Locate the coincidence peak in a 1 ps per bin delay histogram.

  Returns ``(peak_ps, fwhm_ps, floor)``: the maximum of the smoothed histogram,
  its full width at half maximum over the floor, and the floor itself, the
  accidental coincidences per ps.
  """
  smooth = np.convolve(delay_counts, np.ones(smooth_ps) / smooth_ps, mode="same")
  peak = int(np.argmax(smooth))
  floor = float(np.median(smooth))
  above = smooth - floor >= (smooth[peak] - floor) / 2
  left = peak
  while left > 0 and above[left - 1]:
    left -= 1
  right = peak
  while right < above.size - 1 and above[right + 1]:
    right += 1
  return peak, right - left + 1, floor


def select(vault, peak_ps, window_ps=100, bins=512, save=None, progress=None):
  """Second pass: the photons within ``peak_ps +/- window_ps``.

  Returns ``(ghost_image, selected)``. With ``save`` the selected photons are
  also written to that path as a ``.photons`` file.
  """
  count = vault.data.photons.x.shape[0]
  ghost = np.zeros((bins, bins), dtype=np.int64)
  selected = 0
  writer = _Writer(vault, save, peak_ps, window_ps) if save else None
  try:
    for lo, hi in chunks(count):
      delay = delays(vault, lo, hi)
      keep = np.abs(delay - peak_ps) <= window_ps
      x = vault.data.photons.x[lo:hi][keep]
      y = vault.data.photons.y[lo:hi][keep]
      ghost += image_of(x, y, bins)
      if writer:
        writer.append(lo, hi, keep, x, y, delay[keep])
      selected += int(x.size)
      if progress:
        progress(hi, count)
  finally:
    if writer:
      writer.close()
  return ghost, selected


def _dataset(vault, path):
  node = vault.data
  for part in path.strip("/").split("/"):
    node = node[part]
  return node


class _Writer:
  """The selected photons as a ``.photons`` file of the source's own layout.

  Its ``/photons/dt`` is the TDC's delay in ps, one ps a channel, so Preview's
  decay plot draws the coincidence peak as it was resolved.
  """

  COPIED = ("/start/time", "/stop/time")

  def __init__(self, source, path, peak_ps, window_ps):
    if peak_ps + window_ps > np.iinfo(np.uint16).max:
      raise ValueError("a delay past 65535 ps does not fit /photons/dt")
    if os.path.exists(path):
      os.remove(path)
    self._source = source
    self._out = Vault(path, fail_if_exists=True, create_if_missing=True,
                      read_only=False, backend="seven")
    for name, value in source.attributes.items():
      self._out.set_attribute(name, str(value))
    self._out.set_attribute("/photons/TacChannel", "1")
    self._out.set_attribute("/photons/TacOffsetPs", "0")
    self._out.set_attribute("Selection", f"start - stop within {peak_ps} +/- {window_ps} ps")
    self._out.set_attribute("SelectedBy", "ghost.py")
    self._out.set_attribute("SelectedAt", str(datetime.now()))
    self._out.create_dataset("/photons/dt", np.dtype("uint16"))
    for name in ("/photons/x", "/photons/y", "/photons/ms") + self.COPIED:
      self._out.create_dataset(name, np.dtype(_dataset(source, name).dtype))
    # ms[k] is the index of the first photon of millisecond k, so the kept
    # photons need it counted again: how many of them lie ahead of each mark.
    self._ms = source.data.photons.ms[:].astype(np.int64)
    self._new_ms = np.zeros(self._ms.size, dtype=np.uint64)
    self._written = 0

  def append(self, lo, hi, keep, x, y, delay):
    self._out.data.photons.x.append(np.ascontiguousarray(x))
    self._out.data.photons.y.append(np.ascontiguousarray(y))
    self._out.data.photons.dt.append(delay.astype(np.uint16))
    for name in self.COPIED:
      kept = _dataset(self._source, name)[lo:hi][keep]
      _dataset(self._out, name).append(np.ascontiguousarray(kept))
    ahead = np.concatenate(([0], np.cumsum(keep)))
    marks = (self._ms >= lo) & (self._ms < hi)
    self._new_ms[marks] = self._written + ahead[self._ms[marks] - lo]
    self._written += int(keep.sum())
    self._new_ms[self._ms >= hi] = self._written

  def close(self):
    self._out.data.photons.ms.append(self._new_ms)
    self._out.close()


# --- Plots, shared with the notebook ---------------------------------------

def _white(image, percentile):
  """The count drawn white: a percentile, so a few hot pixels do not set it."""
  return max(1.0, float(np.percentile(image, percentile)))


def save_png(path, image, percentile=99.9):
  """Write a count image as an 8-bit grayscale PNG, as ``image_figure`` shows it.

  The same contrast and the same orientation, Preview's: x from right to left.
  """
  shown = np.clip(image[:, ::-1] / _white(image, percentile), 0, 1)
  pixels = np.round(shown * 255).astype(np.uint8)
  height, width = pixels.shape
  # A filter byte ahead of every row, 0 for "none".
  rows = np.hstack((np.zeros((height, 1), dtype=np.uint8), pixels)).tobytes()

  def chunk(kind, data):
    body = kind + data
    return struct.pack(">I", len(data)) + body + struct.pack(">I", zlib.crc32(body))

  with open(path, "wb") as out:
    out.write(b"\x89PNG\r\n\x1a\n")
    out.write(chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 0, 0, 0, 0)))
    out.write(chunk(b"IDAT", zlib.compress(rows, 9)))
    out.write(chunk(b"IEND", b""))

def image_figure(image, title, palette="Greys256", percentile=99.9, size=560):
  """A bokeh figure of a count image, oriented as Photonscore Preview shows it.

  Preview draws y downwards and x from right to left, so both axes run
  backwards here; the counts themselves are left as they are.
  """
  from bokeh.models import ColorBar, LinearColorMapper
  from bokeh.plotting import figure

  mapper = LinearColorMapper(palette=palette, low=0, high=_white(image, percentile))
  fig = figure(title=title, width=size + 70, height=size, match_aspect=True,
               x_range=(POSITION_MAX, 0), y_range=(POSITION_MAX, 0),
               x_axis_label="x", y_axis_label="y",
               tooltips=[("x", "$x{0}"), ("y", "$y{0}"), ("photons", "@image")])
  fig.image(image=[image], x=0, y=0, dw=POSITION_MAX, dh=POSITION_MAX,
            color_mapper=mapper)
  fig.add_layout(ColorBar(color_mapper=mapper, title="photons"), "right")
  fig.grid.visible = False
  return fig


def delay_figure(delay_counts, peak_ps, window_ps, fwhm_ps=None, around_ps=None,
                 width=900, height=360, dark=False):
  """The start-stop histogram, 1 ps a bin, as steps; the selection shaded.

  ``around_ps`` zooms to ``peak_ps +/- around_ps``; without it the whole span
  is drawn. ``dark`` draws it for a dark page: light ink and no background of
  its own, a bokeh plot being unable to follow the page's scheme by itself.
  """
  from bokeh.models import BoxAnnotation, HoverTool, Range1d, Span
  from bokeh.plotting import figure

  lo, hi = 0, delay_counts.size
  if around_ps:
    lo, hi = max(lo, peak_ps - around_ps), min(hi, peak_ps + around_ps)
  # A bin's outline: both of its edges at its own height, so the line steps.
  delay = np.repeat(np.arange(lo, hi + 1), 2)[1:-1]
  photons = np.repeat(delay_counts[lo:hi], 2)
  title = f"Start-stop delay: peak at {peak_ps} ps"
  if fwhm_ps:
    title += f", {fwhm_ps} ps FWHM"
  fig = figure(title=title, width=width, height=height,
               x_axis_label="photon after bucket click, ps",
               y_axis_label="photons per ps")
  fig.y_range = Range1d(0, 1.05 * max(1, int(photons.max())))
  outline = fig.line(delay, photons, line_color="#7fb2e5" if dark else "#4c78a8",
                     line_width=1)
  fig.add_tools(HoverTool(renderers=[outline], mode="vline",
                          tooltips=[("delay", "@x{0} ps"), ("photons", "@y{0}")]))
  fig.add_layout(BoxAnnotation(left=peak_ps - window_ps, right=peak_ps + window_ps,
                               fill_color="#CC3333", fill_alpha=0.3 if dark else 0.18))
  fig.add_layout(Span(location=peak_ps, dimension="height", line_color="#CC3333"))
  if dark:
    _darken(fig)
  return fig


def _darken(fig, ink="#d0d0d0", faint="#444444"):
  """Light ink on no background, for a page that is dark already."""
  fig.background_fill_color = None
  fig.border_fill_color = None
  fig.outline_line_color = faint
  fig.title.text_color = ink
  fig.grid.grid_line_color = faint
  for axis in fig.axis:
    axis.axis_label_text_color = ink
    axis.major_label_text_color = ink
    axis.axis_line_color = ink
    axis.major_tick_line_color = ink
    axis.minor_tick_line_color = ink


def report(path, full, ghost, delay_counts, peak_ps, fwhm_ps, window_ps, title):
  """Write the three plots to one HTML file."""
  from bokeh.io import output_file, save
  from bokeh.layouts import column, row

  output_file(path, title=title)
  save(column(
    row(image_figure(full, f"Every photon: {int(full.sum()):,}"),
        image_figure(ghost, f"Peak +/- {window_ps} ps: {int(ghost.sum()):,}")),
    delay_figure(delay_counts, peak_ps, window_ps, fwhm_ps),
    delay_figure(delay_counts, peak_ps, window_ps, fwhm_ps, around_ps=1000),
  ))


def _progress(label):
  began = time.time()

  def show(done, count):
    print(f"\r  {label}: {100 * done / count:5.1f} %  {time.time() - began:5.1f} s",
          end="\n" if done == count else "", flush=True)
  return show


def main(argv=None):
  parser = argparse.ArgumentParser(
    description="Extract a ghost image from a LINCam .photons file.")
  parser.add_argument("input", help="the .photons file")
  parser.add_argument("--window", type=int, default=100, metavar="PS",
                      help="half width of the selection around the peak (100)")
  parser.add_argument("--bins", type=int, default=RESOLUTIONS[0], choices=RESOLUTIONS,
                      help="pixels along each side of the images (512)")
  parser.add_argument("--save", metavar="FILE",
                      help="also write the selected photons to this .photons file")
  parser.add_argument("--png", action="store_true",
                      help="also write <input>.full.png and <input>.ghost.png")
  parser.add_argument("--html", metavar="FILE",
                      help="where the plots go (<input>.ghost.html)")
  args = parser.parse_args(argv)

  vault = Vault(args.input)
  count = vault.data.photons.x.shape[0]
  print(f"{args.input}: {count:,} photons")

  full, delay_counts = scan(vault, args.bins, progress=_progress("scan"))
  peak_ps, fwhm_ps, floor = find_peak(delay_counts)
  print(f"  peak at {peak_ps} ps, {fwhm_ps} ps FWHM, {floor:.0f} accidentals per ps")

  ghost, selected = select(vault, peak_ps, args.window, args.bins, args.save,
                           progress=_progress("select"))
  accidental = floor * (2 * args.window + 1)
  print(f"  {selected:,} photons within {peak_ps} +/- {args.window} ps "
        f"({100 * selected / count:.2f} %), about {100 * accidental / selected:.1f} % "
        f"of them accidental")
  if args.save:
    print(f"  saved to {args.save}")

  stem = os.path.splitext(args.input)[0]
  html = args.html or stem + ".ghost.html"
  report(html, full, ghost, delay_counts, peak_ps, fwhm_ps, args.window,
         os.path.basename(args.input))
  print(f"  plots in {html}")
  if args.png:
    for name, picture in ((stem + ".full.png", full), (stem + ".ghost.png", ghost)):
      save_png(name, picture)
      print(f"  {name}: {args.bins} x {args.bins}")
  vault.close()


if __name__ == "__main__":
  sys.exit(main())
