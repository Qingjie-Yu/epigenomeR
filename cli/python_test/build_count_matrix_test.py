#!/usr/bin/env python
"""
Test script for build_count_matrix tool that has @tool decorator.
This extracts and calls the underlying function.
"""

from pprint import pprint
import sys
import os

# Add parent directory to path to import from sibling 'tools' folder
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the decorated tool
from python_script.build_count_matrix import build_count_matrix

# Extract the underlying function from the StructuredTool
# The @tool decorator stores the original function in the 'func' attribute
actual_function = build_count_matrix.func

from pathlib import Path

def run_test():
    """Run a test with your BAM files."""

    print("Testing count_matrix_function_with_qc (extracting function from @tool decorator)")
    print("="*60)

    # Example test with BAM files
    bam_dir = Path(
    "/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/data_align/VCOP"
    )
    bam_files = sorted(str(p) for p in bam_dir.glob("*.bam"))
    if not bam_files:
        raise RuntimeError(f"No .bam files found in: {bam_dir}")

    output_dir = "/users/lqi/cli/output"

    print(f"Processing {len(bam_files)} BAM files")
    print(f"Output directory: {output_dir}")
    print("="*60)

    try:
        # Call the extracted function with regions parameter
        result = actual_function(
            bam_paths=bam_files,
            sample="Test",
            output_dir=output_dir,
            regions=800,
        )

        print("\nResult:")
        pprint(result)

        # Check results - build_count_matrix returns {'cm_path': '...'} or {'cm_path': None}
        if result.get('cm_path') is not None:
            print("\nTest PASSED")
            print(f"  Count matrix file: {result['cm_path']}")
            if os.path.exists(result['cm_path']):
                print(f"  File exists: Yes")
                print(f"  File size: {os.path.getsize(result['cm_path'])} bytes")
            else:
                print(f"  WARNING: File does not exist!")
        else:
            print("\nTest FAILED - cm_path is None")

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
    print(f"Type of build_count_matrix: {type(build_count_matrix)}")
    print(f"Has 'func' attribute: {hasattr(build_count_matrix, 'func')}")
    print()

    run_test()
