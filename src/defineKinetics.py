#%%
import os
from src.auxFcts import pickle_load_binary, pickle_dump_binary, print_network_info





def define_kinetics(case_folder):
    system = pickle_load_binary(os.path.join(case_folder, ".NETWORK_system_pickle"))

#%%

