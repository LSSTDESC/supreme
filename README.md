# suprême
**Su**rvey **Pr**op**ê**rty **M**aps with healspars**e** for the LSST survey

This code uses [healsparse](https://github.com/lsstdesc/healsparse) to make
survey property maps for the Vera C. Rubin Observatory Legacy Survey of Space
and Time (LSST) Dark Energy Science Collaboration (DESC).

## Requirements:

`suprême` requires the Rubin LSST stack in order to run, as it is based on the
data butler and stack APIs (including WCS, PSF, and other data).  It is
regularly tested on the stack v19.0.0 and the latest weekly.

In addition to the stack, `suprême` requires
[healsparse](https://github.com/lsstdesc/healsparse) to be installed.

## Installation:

To install the package, go to the parent directory of the source tree and run
`python setup.py install [options]` or use `pip install . [options]`.
Alternatively, you can use EUPS to set up the package with `setup -j -r
/path/to/supreme`.

## Testing:

In order to run the tests, you must also have the
[suprême_testdata](https://github.com/lsstdesc/supreme_testdata) repo.
Assuming you have the stack and `suprême` set up, you can run:

```
setup -j -r /path/to/supreme_testdata
cd tests
pytest test_*rc2.py
pytest test_*dc2.py
```

## Notes:

This software was developed within the LSST DESC using LSST DESC resources, and
so meets the criteria given in, and is bound by, the LSST DESC Publication
Policy for being a “DESC product”.  We welcome requests to access the code for
non-DESC use; if you wish to use the code outside DESC please contact the
developers.

This code is under development and has not yet been released.
