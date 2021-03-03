py_spec (SPEC python utilities)
======

This is the Python package for SPEC.
`py_spec` uses the [black](https://black.readthedocs.io/en/stable/index.html) code style.

Install
======

You can install the package using `pip`.

```
pip install py_spec
```

Or from the GitHub

```
pip install git+https://github.com/PrincetonUniversity/SPEC.git#egg=SPEC&subdirectory=Utilities/pythontools
```

You can also install it locally (`cd /path/to/SPEC/Utlities/pythontools/`).

```
python setup.py install
```
or
```
pip install -e .
```

You can also add the path to your `PYTHONPATH` to directly use this package.
```python
import sys
sys.path.append('/path/to/SPEC/Utlities/pythontools/py_spec')
```

or put this in your env file.
```
export PYTHONPATH=$PYTHONPATH:/path/to/SPEC/Utlities/pythontools/py_spec
```

Use
======

After installation, in your python kernel, you can use the functions by something like
```python
from py_spec import *
```

Documentation
======
You can find documentations at `./docs`.

Developers
======
Here are the contributors:

- Ksenia Aleynikova
- Zhisong Qu
- Jonathan Schilling
- Christopher Berg Smiet
- Caoxiang Zhu