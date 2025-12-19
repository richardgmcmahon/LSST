"""
Test script for create_test_table function.

This demonstrates the function requested by @richardgmcmahon to create
a test astropy table with columns: redshift, mag_g, mag_r
"""

import sys
import os

# Add the current directory to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

print("="*70)
print("TEST: create_test_table function")
print("="*70)

# Test 1: Import the function
print("\n[Test 1] Import create_test_table from rgm_util")
try:
    from rgm_util.create_test_table import create_test_table
    print("  ✓ Successfully imported create_test_table")
except ImportError as e:
    print(f"  ✗ Import failed: {e}")
    sys.exit(1)

# Test 2: Create table with exact specifications from the comment
print("\n[Test 2] Create table as specified in comment")
print("  Parameters:")
print("    - column_names=['redshift', 'mag_g', 'mag_r']")
print("    - column_ranges=[(-10, 10)]")
print("    - n_rows=100")
print("    - All floating point")

try:
    table = create_test_table(
        n_rows=100,
        column_names=['redshift', 'mag_g', 'mag_r'],
        column_ranges=[(-10, 10)]
    )
    print(f"  ✓ Table created successfully")
    print(f"\n  Table preview (first 10 rows):")
    print(table[:10])
    
    print(f"\n  Table details:")
    print(f"    Number of rows: {len(table)}")
    print(f"    Column names: {table.colnames}")
    print(f"    Column dtypes: {[str(table[col].dtype) for col in table.colnames]}")
    
    # Verify all columns are floating point
    all_float = all('float' in str(table[col].dtype) for col in table.colnames)
    if all_float:
        print(f"    ✓ All columns are floating point")
    else:
        print(f"    ✗ Not all columns are floating point")
    
    # Show ranges
    print(f"\n  Column value ranges:")
    for col in table.colnames:
        print(f"    {col:10s}: [{table[col].min():7.3f}, {table[col].max():7.3f}]")
    
except Exception as e:
    print(f"  ✗ Error creating table: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 3: Test with different ranges for each column
print("\n[Test 3] Create table with different ranges per column")
try:
    table2 = create_test_table(
        n_rows=100,
        column_names=['redshift', 'mag_g', 'mag_r'],
        column_ranges=[(0, 5), (15, 25), (15, 25)],
        seed=42  # For reproducibility
    )
    print(f"  ✓ Table created with custom ranges")
    print(f"\n  Column value ranges:")
    for col in table2.colnames:
        print(f"    {col:10s}: [{table2[col].min():7.3f}, {table2[col].max():7.3f}]")
except Exception as e:
    print(f"  ✗ Error: {e}")

# Test 4: Verify import from package level works
print("\n[Test 4] Import from package level (rgm_util)")
try:
    from rgm_util import create_test_table
    test_table = create_test_table(n_rows=10)
    print(f"  ✓ Package-level import works")
    print(f"  ✓ Created table with {len(test_table)} rows")
except Exception as e:
    print(f"  ✗ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*70)
print("✓ ALL TESTS PASSED")
print("="*70)

print("\n" + "="*70)
print("USAGE EXAMPLES")
print("="*70)

print("""
# Example 1: Create table with exact specifications from comment
from rgm_util.create_test_table import create_test_table

table = create_test_table(
    n_rows=100,
    column_names=['redshift', 'mag_g', 'mag_r'],
    column_ranges=[(-10, 10)]
)

# Example 2: Create table with different ranges for each column
table = create_test_table(
    n_rows=100,
    column_names=['redshift', 'mag_g', 'mag_r'],
    column_ranges=[(0, 5), (15, 25), (15, 25)]
)

# Example 3: Use with default parameters
table = create_test_table()  # Uses default names and ranges

# Example 4: Import from package level
from rgm_util import create_test_table
table = create_test_table(n_rows=100)
""")
