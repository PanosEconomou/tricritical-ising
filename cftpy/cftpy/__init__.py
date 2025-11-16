from .cft import *
from .cythonft import * # type: ignore

# Populate the dictionaries with examples listed in SPECIAL_NAMES
for name in SPECIAL_NAMES:
    for LABELS in DICTS:
        for example in SPECIAL_NAMES[name]:
            if name in LABELS:
                LABELS[example] = SPECIAL_NAMES[name][example](LABELS[name])

cft.MODEL_NAMES = list(MODEL_LABELS.keys())

# Assing the appropriate dynamic docstring
for func in DYNAMIC_DOCSTRINGS: fill_docstring(func)
