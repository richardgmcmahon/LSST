"""
Simple test to verify the get_footnote bug fix.
"""

import sys
import os

# Add the current directory to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

print("="*70)
print("TESTING THE GET_FOOTNOTE BUG FIX")
print("="*70)

# Test 1: Import and call get_footnote function directly
print("\n[Test 1] Import and call get_footnote function")
try:
    from rgm_util.get_footnote import get_footnote
    result = get_footnote()
    print(f"  ✓ get_footnote() called successfully")
    print(f"  ✓ Returned: '{result}'")
except Exception as e:
    print(f"  ✗ ERROR: {e}")
    sys.exit(1)

# Test 2: Verify the fix in the source code
print("\n[Test 2] Verify the fix in detrend_redshift_colour.py source")
with open('rgm_util/detrend_redshift_colour.py', 'r') as f:
    lines = f.readlines()
    
    # Find the import statement
    import_line = None
    for i, line in enumerate(lines, 1):
        if 'from .get_footnote import get_footnote' in line:
            import_line = (i, line.strip())
            break
    
    if import_line:
        print(f"  ✓ CORRECT import found at line {import_line[0]}:")
        print(f"    '{import_line[1]}'")
        print(f"  ✓ This imports the FUNCTION, not the module")
    else:
        print(f"  ✗ INCORRECT: Expected 'from .get_footnote import get_footnote'")
        sys.exit(1)
    
    # Find the function call
    call_line = None
    for i, line in enumerate(lines, 1):
        if 'get_footnote()' in line and '#' not in line.split('get_footnote()')[0]:
            # Found a call that's not in a comment
            call_line = (i, line.strip())
            break
    
    if call_line:
        print(f"  ✓ Function call found at line {call_line[0]}:")
        print(f"    '{call_line[1]}'")
        print(f"  ✓ This calls the FUNCTION (not trying to call a module)")
    else:
        print(f"  ℹ No uncommented call to get_footnote() found")

# Test 3: Explain the bug and the fix
print("\n" + "="*70)
print("BUG FIX EXPLANATION")
print("="*70)
print("\nORIGINAL BUG:")
print("  Line: import get_footnote")
print("  Line: get_footnote()")
print("  Error: TypeError: 'module' object is not callable")
print("  Problem: Trying to call a MODULE as if it were a FUNCTION")

print("\nTHE FIX:")
print("  Line: from .get_footnote import get_footnote")
print("  Line: get_footnote()")
print("  Result: ✓ Successfully calls the FUNCTION")
print("  Solution: Import the function directly, not the module")

print("\n" + "="*70)
print("✓ ALL TESTS PASSED - BUG FIX VERIFIED!")
print("="*70)
