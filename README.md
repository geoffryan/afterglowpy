# Semi-analytic GRB Afterglow models

A Python module to interface with Hendrik van Eerten's semi-analytic afterglow model.

## Installation/Building

To build:

```bash
$ python setup.py build_ext -i
```

To install

```bash
$ python setup.py install
```

To grab data from the Open Kilonova Catalog (OKC), ensure `json` and the security setup of `requests` are installed in your local Python environment

```bash
$ pip install requests[security] json
```
