#!/usr/bin/env python
"""
Test script for biclustering_wrapper tool that has @tool decorator.
This extracts and calls the underlying function.
"""

from pprint import pprint
import sys
import os

# Add parent directory to path to import from sibling 'tools' folder
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the decorated tool
from python_script.biclustering_wrapper import biclustering_wrapper

# Extract the underlying function from the StructuredTool
# The @tool decorator stores the original function in the 'func' attribute
actual_function = biclustering_wrapper.func

def run_test():
    """Run a test with your count matrix files."""

    print("Testing biclustering_wrapper")
    print("="*60)

    # Example test with count matrix files
    cm_files = [
        "/users/lqi/cli/test_data/Count_Matrix_800.feather"
    ]

    output_dir = "/users/lqi/cli/output"
    sample = "Test"

    print(f"Processing {len(cm_files)} count matrix file(s)")
    print(f"Output directory: {output_dir}")
    print(f"Sample: {sample}")
    print("="*60)

    try:
        # Call the extracted function
        result = actual_function(
            cm_paths=cm_files,
            sample=sample,
            output_dir=output_dir,
            row_km=15,
            col_km=3
        )

        print("\nResult:")
        pprint(result)

        # Check results - biclustering_wrapper returns {'status': 'success'} or {'status': 'failed'}
        if result.get('status') == 'success':
            print("\nTest PASSED")
            bicluster_dir = os.path.join(output_dir, 'biclustering', sample)
            print(f"  Output directory: {bicluster_dir}")
            if os.path.exists(bicluster_dir):
                print(f"  Directory exists: Yes")
                files = os.listdir(bicluster_dir)
                print(f"  Files created: {len(files)}")
                if files:
                    print(f"  Files: {files}")
            else:
                print(f"  WARNING: Directory does not exist!")
        else:
            print("\nTest FAILED - status is not 'success'")

    except Exception as e:
        print(f"\nTest FAILED with exception: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    # Get parent directory (project root)
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tools_dir = os.path.join(parent_dir, "python_script")

    # Check if we can access the python_script folder
    if not os.path.exists(tools_dir):
        print("ERROR: 'python_script' folder not found!")
        print("Expected location:", tools_dir)
        print("Current script location:", os.path.abspath(__file__))
        sys.exit(1)

    # Check what we're dealing with
    print(f"Type of biclustering_wrapper: {type(biclustering_wrapper)}")
    print(f"Has 'func' attribute: {hasattr(biclustering_wrapper, 'func')}")
    print()

    run_test()
