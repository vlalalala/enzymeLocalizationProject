import os
import re
from pathlib import Path
import signal

class DelayedSignals:
    # sigint is Ctrl+C, snakemake/slurm send SIGTERM and potentially SIGHUP
    SIGNALS_TO_DELAY = [signal.SIGINT, signal.SIGTERM, signal.SIGHUP]

    def __enter__(self):
        self.interrupted = None
        self.old_handlers = {}
        for sig in self.SIGNALS_TO_DELAY:
            self.old_handlers[sig] = signal.signal(sig, self._handler)
        return self

    def _handler(self, sig, frame):
        self.interrupted = sig  # remember which signal came in

    def __exit__(self, type, value, traceback):
        # Restore all original handlers
        for sig, handler in self.old_handlers.items():
            signal.signal(sig, handler)
        # Now raise if we got a signal
        if self.interrupted is not None:
            signal.raise_signal(self.interrupted)

def get_concentrations_files_within_folder(folder):
    """
    Returns a list with the complete paths of files that correspond to
    concentrations (or an empty list, if none exist).
    Compatible with Snakemake wildcards.
    """
    pattern = re.compile(
        r"^interpolation_iteration_nr_\d+_Newton_iteration_nr_\d+_concentrations\.json$"
    )
    files = [
        os.path.join(folder, f)
        for f in os.listdir(folder)
        if pattern.fullmatch(f)
    ]
    return files

def rename_iteration_files(folder, min_digits) -> int:
    files = get_concentrations_files_within_folder(folder)
    folder = Path(folder)
    if not files:
        print(f"No files found in {folder} matching iterations pattern.")
        return min_digits

    files = [Path(f) if os.path.isabs(f) else folder / f for f in files]

    pattern = re.compile(r"(Newton_iteration_nr_)(\d+)")

    # Detect maximum existing digit count in filenames
    max_found_digits = 0
    for f in files:
        match = pattern.search(os.path.basename(f))
        if match:
            max_found_digits = max(max_found_digits, len(match.group(2)))

    final_digits = max(min_digits, max_found_digits)

    for f in files:
        basename = os.path.basename(f)
        match = pattern.search(basename)
        if not match:
            continue

        old_num_str = match.group(2)
        if len(old_num_str) == final_digits:
            continue

        new_num_str = f"{int(old_num_str):0{final_digits}d}"
        new_basename = pattern.sub(lambda m: m.group(1) + new_num_str, basename, count=1)
        old_path = os.path.join(folder, basename)
        new_path = os.path.join(folder, new_basename)

        if old_path != new_path:
            os.rename(old_path, new_path)
            print(f"Renamed: {basename} → {new_basename}")

    return final_digits


