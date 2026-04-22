import os
import shutil
from pathlib import Path
from collections import defaultdict

def copy_and_flatten(root_folder: str, target_filename: str, downloads_dir: str = "Downloads"):
    root_folder = Path(root_folder).resolve()
    downloads_dir = Path(downloads_dir).resolve()

    if not root_folder.exists():
        raise ValueError(f"Root folder does not exist: {root_folder}")

    # Build EXACT mirrored structure inside Downloads
    # e.g. Downloads/data/optimization
    output_base = downloads_dir / Path(*root_folder.parts[-2:])  # <-- FIXED

    name_counts = defaultdict(int)

    for dirpath, _, filenames in os.walk(root_folder):
        for fname in filenames:
            if fname == target_filename:
                full_path = Path(dirpath) / fname

                # relative path strictly INSIDE root_folder
                rel_path = full_path.relative_to(root_folder)

                # flatten path into filename
                flat_name = "_".join(rel_path.parts)

                # deduplicate safely
                name_counts[flat_name] += 1
                count = name_counts[flat_name]

                if count > 1:
                    stem, dot, ext = flat_name.rpartition(".")
                    if dot:
                        flat_name = f"{stem}_{count}.{ext}"
                    else:
                        flat_name = f"{flat_name}_{count}"

                dest_path = output_base / flat_name
                dest_path.parent.mkdir(parents=True, exist_ok=True)

                shutil.copy2(full_path, dest_path)

                print(f"Copied: {full_path} -> {dest_path}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--root_folder", help="e.g. data/optimization")
    parser.add_argument("--filename", help="e.g. interpolation_iteration_nr_0.png")
    parser.add_argument("--downloads", default=".to_download")

    args = parser.parse_args()

    copy_and_flatten(args.root_folder, args.filename, args.downloads)
# python data_search_scripts/search_and_copy.py --root_folder data_private/optimization_enzymaticXtoY_spontaneousYtoZ_1InnerBoundary/ --filename interpolation_iteration_nr_0_final_concentrations.png