```{eval-rst}
.. module:: bindome
```

```{eval-rst}
.. automodule:: bindome
    :noindex:
```

# API

```
import bindome as bd
```

Once the module is imported, the retrieval of datasets is based on the reference sources, and protein names

E.g.

```
bd.constants.ANNOTATIONS_DIRECTORY = "/your/path/to/annotations"
peaks = mb.bindome.datasets.REMAP2020.get_remap_peaks("GATA1")
```

Dataset examples:

- ChIP-atlas
- ReMap
- SELEX
- PBM
- ProBound
- scATAC

## Tools: `tl`

```{eval-rst}
.. module:: bindome.tl
.. currentmodule:: bindome

.. autosummary::
    :toctree: generated/

```

```{eval-rst}
.. autosummary::
    :toctree: generated/
```

## Plotting

```{eval-rst}
.. module:: bindome.pl
.. currentmodule:: bindome

.. autosummary::
    :toctree: generated/
```
