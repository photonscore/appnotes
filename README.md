# Photonscore application notes

The code behind the application notes on
[photonscore.de](https://www.photonscore.de): for each note, the script that
reads the measured data and makes the note's figures, and a notebook that walks
through it step by step. The data files themselves are linked from each note;
they run to gigabytes and are not kept here.

| | note | instrument | what the code does |
| --- | --- | --- | --- |
| [1](1/) | [Quantum ghost imaging](https://www.photonscore.de/appnotes/1) | LINCam | picks the photon pairs out of 244 million single photons by their picosecond time stamps, and draws the sample from the ones that never met it |

## How a note is laid out

A folder per note, named by the note's number on the site:

```
1/
  README.md      what was measured, how to install, where the data is, how to run
  ghost.py       the script: a command line tool and the functions the notebook calls
  ghost.ipynb    the same steps one at a time, a plot after each
  result.png     what comes out
```

Every note's README has its own install line. They all need
[uv](https://docs.astral.sh/uv/getting-started/installation/), which fetches
Python as well, and the `photonscore` Python package, which reads and writes
`.photons` files and comes from the
[downloads page](https://www.photonscore.de/support/downloads). Files a script
writes open in Photonscore Preview, from the same page.

Questions about a note or an instrument:
[email@photonscore.de](mailto:email@photonscore.de).

## License

[MIT](LICENSE).
