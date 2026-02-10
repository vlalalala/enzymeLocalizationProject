"""
To create .csv templates, run from snakemake root directory
src/create_file_structure.py

"""
#%%
from __future__ import annotations

import sys
from pathlib import Path
from itertools import product
import shutil


def build_options_cartesian_product(
    root_dir: str | Path,
    output_dir: str | Path,
    options_prefix: str = "options_",
    combo_prefix: str = "combo_",
    max_folders: int | None = None,
) -> int:
    """
    Create folders with all combinations built from combo_* files across options_* folders.

    Each output folder contains exactly one file from each options_* folder.
    Copied files are renamed:
        basename = suffix of options_* folder (the '*')
        extension = original extension of the chosen combo file

    Example:
        options_alpha/combo_000001.csv -> alpha.csv
        options_beta/combo_000010.yaml -> beta.yaml

    Returns number of output folders created.
    Function written by ChatGPT.
    """
    root = Path(root_dir)
    out_root = Path(output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    # Find options_* folders
    options_dirs = sorted(
        [p for p in root.iterdir() if p.is_dir() and p.name.startswith(options_prefix)],
        key=lambda p: p.name,
    )
    if not options_dirs:
        raise ValueError(f"No '{options_prefix}*' folders found under: {root}")

    # For each options_* folder, gather combo_* files (any extension)
    bundles: list[tuple[str, list[Path]]] = []
    for opt_dir in options_dirs:
        suffix = opt_dir.name[len(options_prefix):]
        if suffix == "":
            raise ValueError(f"Options folder has empty suffix: {opt_dir.name}")

        combo_files = sorted(
            [f for f in opt_dir.iterdir()
             if f.is_file() and f.name.startswith(combo_prefix)],
            key=lambda p: p.name,
        )
        if not combo_files:
            raise ValueError(f"No '{combo_prefix}*' files found in: {opt_dir}")

        bundles.append((suffix, combo_files))

    # Cartesian product over file lists
    lists_only = [files for (_suffix, files) in bundles]

    created = 0
    for idx, selection in enumerate(product(*lists_only), start=1):
        if max_folders is not None and created >= max_folders:
            break

        folder_name = f"combined_{idx:06d}"
        combo_out_dir = out_root / folder_name
        combo_out_dir.mkdir(parents=False, exist_ok=False)

        # Copy each selected file, renaming to "<options_suffix><original_extension>"
        for (suffix, _files), chosen in zip(bundles, selection):
            ext = chosen.suffix  # includes leading dot, e.g. ".csv" / ".yaml"
            new_name = f"{suffix}{ext}"
            dest = combo_out_dir / new_name
            shutil.copy2(chosen, dest)

        created += 1

    return created


if __name__ == "__main__":
    folder_with_parameter_ranges = sys.argv[1]
    n = build_options_cartesian_product(
        root_dir=folder_with_parameter_ranges,              # where options_* folders live
        output_dir=folder_with_parameter_ranges,    # where combined_* folders will be created
        max_folders=None,           # set to an int to cap output
    )
    print(f"Created {n} combined folders.")