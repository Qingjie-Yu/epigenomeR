#!/usr/bin/env python
"""
Test script for qc_by_percentile tool that has @tool decorator.
This extracts and calls the underlying function.
"""

from pprint import pprint
import sys
import os

# Add parent directory to path to import from sibling 'tools' folder
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the decorated tool
from python_script.qc_by_percentile import qc_by_percentile

# Extract the underlying function from the StructuredTool
# The @tool decorator stores the original function in the 'func' attribute
actual_function = qc_by_percentile.func

from pathlib import Path

def run_test():
    """Run a test with BAM/BED files."""

    print("Testing qc_by_percentile (extracting function from @tool decorator)")
    print("="*60)

    # Example test with BAM/BED files
    bam_dir = Path(
    "/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/data_align/VCOP"
    )
    input_files = sorted(str(p) for p in bam_dir.glob("*.bam"))
    if not input_files:
        raise RuntimeError(f"No .bam files found in: {bam_dir}")

    output_dir = "/users/lqi/cli/output"
    sample = "Test"

    print(f"Processing {len(input_files)} file(s)")
    print(f"Output directory: {output_dir}")
    print("="*60)

    try:
        # Call the extracted function
        result = actual_function(
            file_paths=input_files,
            sample=sample,
            output_dir=output_dir,
            filtered_percentile=0.25
        )

        print("\nResult:")
        pprint(result)

        # Check results - qc_by_percentile returns {'status': 'success'} or {'status': 'failed'}
        if result.get('status') == 'success':
            print("\nTest PASSED")
            qc_dir = os.path.join(output_dir, 'qc', sample)
            print(f"  Output directory: {qc_dir}")
            if os.path.exists(qc_dir):
                print(f"  Directory exists: Yes")
                files = os.listdir(qc_dir)
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
    print(f"Type of qc_by_percentile: {type(qc_by_percentile)}")
    print(f"Has 'func' attribute: {hasattr(qc_by_percentile, 'func')}")
    print()

    run_test()
