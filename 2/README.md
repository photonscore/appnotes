# 2 · Phasor analysis of LINCam data with PhasorPy

![The photons of a Convallaria section, their phasor plot, the pixels of each
cluster on it, and the share of the shorter of two lifetimes](result.png)

*Top left: the 284 million photons of a* Convallaria majalis *section,
recorded with a LINCam in 23 minutes. Top right: their phasors at 59.5 MHz,
one point a pixel, with three clusters found on the plot and the line between
the two lifetimes, 0.64 and 2.92 ns, whose mixtures the cloud follows. Bottom
left: the pixels of each cluster, in its colour. Bottom right: the share of
the 0.64 ns component in every pixel.*

A LINCam records every photon with its position and its arrival after the
laser pulse. The `photonscore` package reads the `.photons` file, and from the
first histogram on everything is [PhasorPy](https://www.phasorpy.org), an
open-source library for phasor analysis: no fitting and no model chosen in
advance, only the Fourier coefficients of each pixel's decay and what the
phasor plot makes of them. The application note itself:
[photonscore.de/appnotes/2](https://www.photonscore.de/appnotes/2).

## Install

Get [uv](https://docs.astral.sh/uv/getting-started/installation/), then paste
one line in this folder. It fetches Python 3.12 as well, PhasorPy, and the
`photonscore` package as
[photonscore.de](https://www.photonscore.de/support/downloads) has it.

Windows:

```powershell
uv venv --python 3.12; uv pip install jupyterlab matplotlib phasorpy==0.12 https://www.photonscore.de/support/downloads/26.0.72/photonscore-26.0.72-windows-py3.12.zip
```

macOS:

```bash
uv venv --python 3.12; uv pip install jupyterlab matplotlib phasorpy==0.12 https://www.photonscore.de/support/downloads/26.0.72/photonscore-26.0.72-mac-py3.12.zip
```

Linux:

```bash
uv venv --python 3.12; uv pip install jupyterlab matplotlib phasorpy==0.12 https://www.photonscore.de/support/downloads/26.0.72/photonscore-26.0.72-linux-py3.12.zip
```

## Data

- [convallaria.photons](https://www.photonscore.de/appnotes/2/data/convallaria.photons)
  — the recording: 284,223,008 photons over 1371 s, 1.6 GB.

## Run

```bash
uv run phasor.py convallaria.photons
```

Some 15 seconds on a laptop: it streams the file in chunks, so the memory in
use stays near 2 GB whatever the file's length. It prints what it found and
writes `convallaria.phasor.png`, the figure above, and
`convallaria.phasor.ome.tif`, the calibrated phasors of the first three
harmonics.

| option | |
| --- | --- |
| `--pixels N` | pixels along each side of the image, 512 unless said |
| `--bins N` | time bins over the 4096 TAC channels, 256 unless said |
| `--harmonic N` | the harmonic plotted and analysed: 1, 2 or 3, 3 unless said |
| `--min-photons N` | pixels with fewer photons are left out, 300 unless said |
| `--clusters N` | populations sought on the phasor plot, 3 unless said |
| `--dark` | draw the figure for a dark page, as `<input>.phasor-dark.png` |

The notebook walks through the same steps, a plot after each, and shows the
phasors before calibration and at all three harmonics:

```bash
uv run jupyter lab phasor.ipynb
```

## What the script does

1. **Histograms the photons** into a (y, x, time) cube, 512 × 512 pixels and
   256 time bins, and the decay of all photons at the file's own 12.3 ps a
   channel.
2. **Takes the phasor of every pixel** with PhasorPy's `phasor_from_signal`.
   The 4096 channels span 50.4 ns, and the decay dies away well inside them,
   so the histogram is one period of a periodic signal and its harmonics are
   19.8, 39.7 and 59.5 MHz.
3. **Calibrates** with `phasor_calibrate`. The reference is the instrument
   response: some 80 ps wide, far narrower than a period, its phasor is a
   point on the unit circle at the phase of time zero, a reference of
   lifetime 0. Time zero is where the decay of all photons climbs fastest,
   10.523 ns in this recording; on simulated photons this lands within 3 ps
   of the truth, and a known single exponential comes out within 0.0005 of
   its place on the semicircle.
4. **Leaves out the dim pixels** with `phasor_threshold` and median-filters
   the rest with `phasor_filter_median`: 74,426 pixels have 300 photons or
   more.
5. **Finds three populations** on the phasor plot with `phasor_cluster_gmm`,
   whose ellipses serve as cursors: `mask_from_elliptic_cursor` and
   `pseudo_color` colour their pixels. In phase and modulation lifetime at
   59.5 MHz they are 0.99 and 1.78 ns, 1.52 and 2.22 ns, and 1.93 and 2.68 ns.
   All three lie inside the semicircle: no pixel decays with a single
   lifetime.
6. **Splits every pixel into two lifetimes.** The cloud's long axis meets the
   semicircle, `phasor_semicircle_intersect`, at 0.64 and 2.92 ns, and
   `phasor_component_fraction` places each pixel along the line between them.
   The cloud is wider than a line, so more than two lifetimes are at work and
   the share is a summary of each pixel, not a fit.

What it reads from the file:

| dataset | | |
| --- | --- | --- |
| `/photons/x`, `/photons/y` | `uint16` | position, 12 bit |
| `/photons/dt` | `uint16` | arrival after the laser pulse, 12 bit |
| attribute `/photons/TacChannel` | | ps a `dt` channel, 12.3 here |

Images are drawn the way Photonscore Preview draws them, x from right to left.
