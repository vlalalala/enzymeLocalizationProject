import argparse
from pathlib import Path
from auxiliary_functions import dump_json
from auxiliary_functions_using_standard_library import load_json


def create_convergence_file():
    # Create file with information on convergence
    subfolders = sorted([
    p for p in Path(FOLDER_TO_STUDY).iterdir()
    if p.is_dir() and p.name.startswith("combined_")
    ])

    convergence_data = {}
    convergence_prefix = "with early convergence: "

    for subfolder in subfolders:
        key = subfolder.name.split("combined_", 1)[1]  # everything after "combined_"
        log_file = subfolder / ".newton_solver.log"
        concentrations_file = subfolder / ".species_steady_state_concentrations.json"

        if log_file.exists():
            # Read last line from log
            with open(log_file, "r") as f:
                lines = f.readlines()
                last_line = lines[-1].strip() if lines else ""
            # Remove prefix if present
            if last_line.startswith(convergence_prefix):
                last_line = last_line[len(convergence_prefix):]
                value = last_line
            else:
                value = "running simulation (not finished yet)"

        elif concentrations_file.exists():
            # this case should only occur if the simulation should not run
            concentrations_dict = load_json(concentrations_file)
            value = f"simulation not run: {concentrations_dict["error"]}"
        else:
            value = "neither log nor concentrations file found"
        convergence_data[key] = value
    dump_json(FOLDER_TO_STUDY, "convergence_data", convergence_data)

if __name__ == "__main__":
    # Parse arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_to_study", type=str, help="Path to folder with system info")
    args = parser.parse_args()

    FOLDER_TO_STUDY = args.folder_to_study

    create_convergence_file()
    

