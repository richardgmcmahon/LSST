#!/usr/bin/env python3
"""
Example script demonstrating create_test_table function.

This creates a test astropy table with the exact specifications requested:
- column_names=['redshift', 'mag_g', 'mag_r']
- column_ranges=[(-10, 10)]
- n_rows=100
- All floating point
"""

from rgm_util.create_test_table import create_test_table

# Create the table with exact specifications from the request
table = create_test_table(
    n_rows=100,
    column_names=['redshift', 'mag_g', 'mag_r'],
    column_ranges=[(-10, 10)]
)

print("="*70)
print("Test Astropy Table Created")
print("="*70)
print(f"\nTable Info:")
print(f"  Rows: {len(table)}")
print(f"  Columns: {table.colnames}")
print(f"  Data types: {[str(table[col].dtype) for col in table.colnames]}")

print(f"\nFirst 10 rows:")
print(table[:10])

print(f"\nColumn Statistics:")
for col in table.colnames:
    print(f"  {col}:")
    print(f"    Min:  {table[col].min():8.3f}")
    print(f"    Max:  {table[col].max():8.3f}")
    print(f"    Mean: {table[col].mean():8.3f}")
    print(f"    Std:  {table[col].std():8.3f}")

print("\n" + "="*70)
print("Table is ready to use!")
print("="*70)
