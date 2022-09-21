# bindome

[![Tests][badge-tests]][link-tests]
[![Documentation][badge-docs]][link-docs]

[badge-tests]: https://img.shields.io/github/workflow/status/ilibarra/bindome/Test/main
[link-tests]: https://github.com/theislab/bindome/actions/workflows/test.yml
[badge-docs]: https://img.shields.io/readthedocs/bindome

## Getting started

Please refer to the [documentation][link-docs]. In particular, the

-   [API documentation][link-api].

## What is bindome ?

Bindome is a repository that gathers and assembles biomolecular binding data (protein-DNA/RNA) from genomics studies into data representations that allow downstream ML-tasks.

### How to use it?

Import the package

```python
import bindome as bd
```

## Installation

You need to have Python 3.8 or newer installed on your system. If you don't have
Python installed, we recommend installing `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`\_.

There are several alternative options to install bindome:

<!--
1) Install the latest release of `bindome` from `PyPI <https://pypi.org/project/bindome/>`_:

```bash
pip install bindome
```
-->

1. Install the latest development version:

```bash
pip install git+https://github.com/theislab/bindome.git@main
```

## Release notes

See the [changelog][changelog].

## Contact

For questions and help requests, you can reach out in the [scverse discourse][scverse-discourse].
If you found a bug, please use the [issue tracker][issue-tracker].

## Citation

for cite bindome, please use the following:

```bibtex
@software{bindome,
  author = {Ibarra}},
  doi = {},
  month = {},
  title = {{bindome}},
  url = {https://github.com/theislab/bindome},
  year = {2022}
}
```

[scverse-discourse]: https://discourse.scverse.org/
[issue-tracker]: https://github.com/theislab/bindome/issues
[changelog]: https://bindome.readthedocs.io/en/latest/changelog.html
[link-docs]: https://bindome.readthedocs.io/en/latest/#
[link-api]: https://bindome.readthedocs.io/en/latest/api.html
