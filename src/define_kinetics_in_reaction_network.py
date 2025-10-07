#%%
import pickle
from check_validity_chemical_network import check_validity_of_csv_files


def pickle_load_binary(path):
    with open(path, 'rb') as f:
        loaded_variable = pickle.load(f)
    return loaded_variable
#%%
# Get the absolute path to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Go up to the project root (one or more levels as needed)
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))


violacein_folder = os.path.join(
    PROJECT_ROOT, *["data", "violacein_0"]
)

#%%
def define_differential_equations(case_folder):
    system = pickle_load_binary(os.path.join(case_folder, ))

#%%
define_differential_equations(violacein_folder)

# %%
