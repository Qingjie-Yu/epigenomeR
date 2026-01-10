import os
import subprocess
from typing import Dict, Optional, List, Union
from langchain.tools import tool
from .LogCapture import auto_log_capture

@tool
@auto_log_capture()
def biclustering_wrapper(
    cm_paths: Union[str, List[str]],
    sample: str,
    output_dir: str,
    no_filter: bool = False,
    no_annotation: bool = False,
    no_plot: bool = False,
    ref_genome: str = "hg38",
    ref_source: str = "knownGene",
    distributions: Union[str, List[str]] = "genic,ccre",
    row_km: int = 15,
    col_km: int = 3
) -> Dict:
    """Perform biclustering analysis on chromatin accessibility count matrices.

    This tool performs biclustering analysis on count matrices, including optional filtering,
    k-means clustering, and gene annotation of genomic regions. The process involves:
    1. Loading count matrix files in Feather format
    2. Optional filtering for highly variable regions (HVR)
    3. K-means clustering on rows (regions) and columns (samples)
    4. Optional gene annotation and TFBS analysis
    5. Generating visualization plots

    Args:
        cm_paths: Absolute path(s) to input count matrix files (.feather). Can be:
                 - A single path string: "/path/to/matrix.feather"
                 - A list of paths: ["/path/to/matrix1.feather", "/path/to/matrix2.feather"]
        sample: Name of the sample being processed. Used for creating output subdirectory.
        output_dir: Base output directory where results will be saved.
        no_filter: If True, disable HVR filtering. Default: False.
        no_annotation: If True, disable annotation and TFBS steps. Default: False.
        no_plot: If True, disable plotting. Default: False.
        ref_genome: Reference genome assembly. Options: "hg38" or "mm10". Default: "hg38".
        ref_source: Annotation source. Options: "knownGene" or "GENCODE". Default: "knownGene".
        distributions: Comma-separated distribution types or list. Default: "genic,ccre".
        row_km: Number of k-means clusters for rows. Default: 15.
        col_km: Number of k-means clusters for columns. Default: 3.

    Returns:
        Dict: Dictionary containing:
            - 'status': 'success' if completed successfully, 'failed' if error occurred.

    Example:
        result = biclustering_wrapper(
            cm_paths="/output/count_matrix/C1/C1_Count_Matrix_800.feather",
            sample="C1",
            output_dir="/output",
            row_km=15,
            col_km=3
        )
        # Returns: {'status': 'success'}
    """

    # Create output directory with sample subdirectory
    bicluster_dict = os.path.join(output_dir, 'biclustering', sample)
    os.makedirs(bicluster_dict, exist_ok=True)

    # Get R script path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path = os.path.join(current_dir, "R_scripts", "biclustering_wrapper_cli.R")

    # Convert cm_paths to comma-separated string
    cm_paths_str = ",".join(cm_paths) if isinstance(cm_paths, list) else str(cm_paths)

    # Convert distributions to string if list
    distributions_str = ",".join(distributions) if isinstance(distributions, list) else str(distributions)

    # Build command
    cmd = [
        "Rscript", r_script_path,
        "--cm_paths", cm_paths_str,
        "--out_dir", bicluster_dict,
        "--ref_genome", ref_genome,
        "--ref_source", ref_source,
        "--distributions", distributions_str,
        "--row_km", str(row_km),
        "--col_km", str(col_km)
    ]

    if no_filter:
        cmd.append("--no_filter")
    if no_annotation:
        cmd.append("--no_annotation")
    if no_plot:
        cmd.append("--no_plot")

    try:
        subprocess.run(cmd, check=True, cwd=os.path.dirname(r_script_path))
        return {"status": "success"}
    except Exception as e:
        print(e)
        return {"status": "failed"}
