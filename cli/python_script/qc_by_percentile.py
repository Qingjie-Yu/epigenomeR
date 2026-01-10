import os
import subprocess
from typing import Dict, Optional, List, Union
from langchain.tools import tool
from .LogCapture import auto_log_capture

@tool
@auto_log_capture()
def qc_by_percentile(
    file_paths: Union[str, List[str]],
    sample: str,
    output_dir: str,
    filtered_percentile: float = 0.25,
    no_plot: bool = False,
    split_pair_by: str = "-",
    group_csv: Optional[str] = None,
    crf_col: str = "crf",
    category_col: str = "category"
) -> Dict:
    """Perform QC on BAM/BED files using percentile filtering.

    This tool performs quality control on BAM or BED files by filtering based on
    percentile thresholds and optionally generates visualization heatmaps. The process:
    1. Reads BAM/BED files
    2. Applies percentile-based filtering to identify high-quality regions
    3. Optionally generates heatmap visualizations
    4. Saves filtered results and QC metrics

    Args:
        file_paths: Absolute path(s) to input BAM or BED files. Can be:
                   - A single path string: "/path/to/file.bam"
                   - A list of paths: ["/path/to/file1.bam", "/path/to/file2.bed"]
        sample: Name of the sample being processed. Used for creating output subdirectory.
        output_dir: Base output directory where QC results will be saved.
        filtered_percentile: Percentile threshold for QC filtering. Default: 0.25.
        no_plot: If True, disable heatmap generation. Default: False.
        split_pair_by: Delimiter for splitting CRF pairs. Default: "-".
        group_csv: Optional CSV file with CRF grouping information. Default: None.
        crf_col: Column name for CRF identifiers in group CSV. Default: "crf".
        category_col: Column name for CRF categories in group CSV. Default: "category".

    Returns:
        Dict: Dictionary containing:
            - 'status': 'success' if completed successfully, 'failed' if error occurred.

    Example:
        result = qc_by_percentile(
            file_paths=["/data/sample1.bam", "/data/sample2.bam"],
            sample="H3K27ac",
            output_dir="/output",
            filtered_percentile=0.25
        )
        # Returns: {'status': 'success'}
    """

    # Create output directory with sample subdirectory
    qc_dir = os.path.join(output_dir, 'qc', sample)
    os.makedirs(qc_dir, exist_ok=True)

    # Get R script path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path = os.path.join(current_dir, "R_scripts", "qc_by_percentile_cli.R")

    # Convert file_paths to comma-separated string
    file_paths_str = ",".join(file_paths) if isinstance(file_paths, list) else str(file_paths)

    # Build command
    cmd = [
        "Rscript", r_script_path,
        "--file_paths", file_paths_str,
        "--out_dir", qc_dir,
        "--filtered_percentile", str(filtered_percentile),
        "--split_pair_by", split_pair_by,
        "--crf_col", crf_col,
        "--category_col", category_col
    ]

    if no_plot:
        cmd.append("--no_plot")
    if group_csv:
        cmd.extend(["--group_csv", group_csv])

    try:
        subprocess.run(cmd, check=True, cwd=os.path.dirname(r_script_path))
        return {"status": "success"}
    except Exception as e:
        print(e)
        return {"status": "failed"}
