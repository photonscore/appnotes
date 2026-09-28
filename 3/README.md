# 3 · Lifetime fitting of LINCam data with FLIMKit

![The photons of a Convallaria section, their decay with a four-lifetime
fit, the mean lifetime of every pixel, and how the photons spread over
it](result.png)

*Top left: the 284 million photons of the* Convallaria majalis *section of
note 2. Top right: their decay, 49 ps a bin, with FLIMKit's reconvolution fit
of four lifetimes, 0.20, 0.81, 2.55 and 5.12 ns, and the residuals below.
Bottom: every pixel fitted with those four lifetimes held, 2 × 2 pixels
binned, as its intensity-weighted mean lifetime, and how the photons spread
over it: three populations, at 1.4, 2.0 and 2.6 ns.*

[FLIMKit](https://github.com/FLIMKit/FLIMKit) is an open-source toolkit for
fluorescence lifetime imaging: reconvolution fitting, lifetime distributions,
phasors, a desktop application and a Python API. It opens LINCam `.photons`
files itself, through [photonsfile](https://github.com/alex1075/photonsfile),
a pure-Python reader written from the
[D7 format](https://github.com/photonscore/d7) Photonscore opened. This note
fits the recording of [note 2](../2/) with it, and the two agree: in the three
populations note 2 found on its phasor plot, FLIMKit's mean lifetime is 1.45,
2.01 and 2.59 ns, each between the phase and the modulation lifetime PhasorPy
gives the same pixels. The application note itself:
[photonscore.de/appnotes/3](https://www.photonscore.de/appnotes/3).

## Install

Get [uv](https://docs.astral.sh/uv/getting-started/installation/), then paste
one line in this folder, on Windows, macOS or Linux alike. It fetches Python
3.12 as well, and FLIMKit with photonsfile.

```bash
uv venv --python 3.12; uv pip install jupyterlab flimkit==0.13.6
```

## Data

- [convallaria.photons](https://www.photonscore.de/appnotes/2/data/convallaria.photons)
  — the recording, note 2's download: 284,223,008 photons over 1371 s,
  1.6 GB.
- [convallaria-40M.photons](https://www.photonscore.de/appnotes/3/data/convallaria-40M.photons)
  — its first 225 s, 40,000,500 photons, 222 MB: small enough for FLIMKit to
  open in some 3.5 GB of memory, on a laptop or in its desktop application.

## Run

```bash
uv run fit.py convallaria.photons
```

FLIMKit reads the whole file into memory, some 80 bytes a photon: the 284
million photons take 24 GB, so a machine with 32 GB. On a smaller one, fit the
first 225 s instead, with the pixel threshold lowered to match:

```bash
uv run fit.py convallaria-40M.photons --min-photons 150
```

It prints the fit and
writes `convallaria.fit.png`, the figure above, and `convallaria.fit.tif`,
the maps as 32-bit float planes: intensity, the two mean lifetimes, the
reduced χ², and each lifetime's amplitude and share of the photons.

| option | |
| --- | --- |
| `--nexp N` | exponentials in the model, 1 to 4, 4 unless said |
| `--fit-start NS` | where the fit starts, half the peak unless said |
| `--irf-fwhm NS` | width of the Gaussian response, 0.08 ns unless said |
| `--bins N` | time bins over the 4096 TAC channels, 1024 unless said |
| `--pixels N` | pixels along each side of the image, 512 unless said |
| `--binning N` | pixels summed along each side for the fit, 2 unless said |
| `--min-photons N` | pixels with fewer photons are left out, 1000 unless said |
| `--dark` | draw the figure for a dark page, as `<input>.fit-dark.png` |

The notebook walks through the same steps, a plot after each, and shows what
one to four exponentials make of the decay:

```bash
uv run jupyter lab fit.ipynb
```

FLIMKit's desktop application opens the same files; its form takes up to
three exponentials, the Python API four. FLIMKit's own documentation is on its
[wiki](https://github.com/FLIMKit/FLIMKit/wiki).

## What the script does

1. **Reads the file with FLIMKit**, `FLIMFile`: every photon binned by its
   position, 512 × 512 pixels, and by its arrival after the laser pulse,
   1024 bins of 49.2 ps over the 4096 TAC channels. The time axis comes from
   the file's `/photons/TacChannel`, 12.3 ps a channel.
2. **Models the instrument response** as a Gaussian of 80 ps FWHM at the
   decay's peak, `gaussian_irf_from_fwhm`; the fit shifts it and may widen
   it. LINCam's response has a core some 70 ps wide and, ahead of it, a weak
   shoulder, a thousandth of the peak and 1 ns long, which no Gaussian
   follows. So the fit starts where the decay climbs through half its peak,
   10.53 ns in this recording.
3. **Fits the decay of all photons**, `fit_summed`: four exponentials
   convolved with the response, plus a flat background, on the Poisson
   likelihood, by differential evolution and a Levenberg-Marquardt polish.
   The lifetimes are 0.20, 0.81, 2.55 and 5.12 ns, in amplitude 46, 32, 20
   and 1.5 %, and the intensity-weighted mean is 2.06 ns. Each exponential
   added cuts the reduced χ² about tenfold, 33,700, 1,810 and 203 for one to
   three, down to 15.3 for four, while the mean lifetime stays at 2.05 to
   2.06 ns.
4. **Fits every pixel**, `fit_per_pixel`, with the four lifetimes held and
   only their amplitudes free, solved by non-negative least squares over the
   same window: 256 × 256 pixels after 2 × 2 binning, of which 19,653 have
   1000 photons or more and hold 98 % of the photons.
5. **Maps the intensity-weighted mean lifetime** and histograms it, each
   pixel weighted by its photons.

The residuals show one feature the fit leaves: a small peak 5.9 ns after the
first, at 16.4 ns. FLIMKit can leave such a band out of the fit with
`exclude_ns`; it is kept in here.

What FLIMKit reads from the file:

| dataset | | |
| --- | --- | --- |
| `/photons/x`, `/photons/y` | `uint16` | position, 12 bit |
| `/photons/dt` | `uint16` | arrival after the laser pulse, 12 bit |
| attribute `/photons/TacChannel` | | ps a `dt` channel, 12.3 here |

Images are drawn the way Photonscore Preview draws them, x from right to left.
