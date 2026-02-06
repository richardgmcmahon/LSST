# LSST Data Processing Utilities

This repository contains utilities and examples for working with LSST data and astronomical catalogs.

## VizieR Table Dimensions

### How to get the number of rows and columns in a VizieR table using astroquery

When using `astroquery v0.4.12.dev255` with the `VizierClass`, you can easily get the dimensions of a returned table:

#### Quick Answer

```python
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u

# Create VizieR query
v = Vizier()
v.ROW_LIMIT = 100

# Query a catalog
coord = SkyCoord(ra=10.68, dec=41.27, unit='deg')
result = v.query_region(coord, radius=5*u.arcmin, catalog="II/246")

if result:
    table = result[0]
    n_rows = len(table)              # Number of rows
    n_cols = len(table.colnames)     # Number of columns
    print(f"Rows: {n_rows}, Columns: {n_cols}")
```

#### Using Helper Functions

This repository provides `vizier_utils.py` with convenient helper functions:

```python
from astroquery.vizier import Vizier
from vizier_utils import get_vizier_table_dimensions, print_vizier_table_info

v = Vizier()
result = v.query_object("M31", catalog=["II/246"])

if result:
    table = result[0]
    
    # Get dimensions
    n_rows, n_cols = get_vizier_table_dimensions(table)
    print(f"Rows: {n_rows}, Columns: {n_cols}")
    
    # Or print comprehensive info
    print_vizier_table_info(table, name="2MASS M31")
```

#### Key Methods

- `len(table)` - Returns the number of rows
- `len(table.colnames)` - Returns the number of columns
- `table.colnames` - Returns a list of all column names
- `table.info()` - Shows detailed table information including column types

#### Complete Example

See `vizier_table_dimensions_example.py` for a complete working example with multiple query methods:

```bash
python vizier_table_dimensions_example.py
```

The example demonstrates:
- Basic queries to get table dimensions
- Querying with specific columns
- Querying by object name
- Handling multiple tables in results

#### Installation

To use these examples, you need astroquery:

```bash
pip install astroquery
```

For the specific version mentioned:
```bash
pip install astroquery==0.4.12.dev255
```

Or for the latest stable version:
```bash
pip install --upgrade astroquery
```

## Repository Contents

- `vizier_table_dimensions_example.py` - Examples for getting VizieR table dimensions
- `vizier_utils.py` - Utility functions for working with VizieR tables
- `explore_rsp_tap.py` - TAP service exploration utilities
- `query_dp1.py` - DP1 data queries
- `lsst_util.py` - LSST utility functions
- `xmatch_dp1.py` - Cross-matching utilities

## Additional Resources

- [astroquery VizieR documentation](https://astroquery.readthedocs.io/en/latest/vizier/vizier.html)
- [VizieR catalog access](http://vizier.u-strasbg.fr/)
- [LSST DP1 documentation](https://sdm-schemas.lsst.io/dp1.html)
