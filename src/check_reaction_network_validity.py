#%%
import sys
import os
import pandas as pd
from auxiliary_functions_using_standard_library import load_json, pickle_dump_binary
from auxiliary_functions import (
    define_ratio_from_string, no_empty_cells, no_repeated_rows_in_csv_file,
    checks_lack_of_repetitions, check_correct_type, define_region_list)

def create_pandas_dataframe_from_csv_file(csv_file: str, num_regions: int) -> pd.DataFrame:
    """Reads the csv files and creates corresponding panda dataframes.
    In case some columns have "ratio" as a substring of the header, the values
    are converted to Ratio types (or an error is raised if that is not possible).
    The types for each entry in each column are checked.
    """
    assert no_repeated_rows_in_csv_file(csv_file)
    dataframe = pd.read_csv(csv_file, encoding='utf-8-sig')
    # Checks species and enzymes
    string_cols = dataframe.columns[dataframe.columns.str.contains(r"species|enzyme|name")]
    assert check_correct_type(dataframe, string_cols, str)
    # Checks k, hill, constants
    float_cols = dataframe.columns[dataframe.columns.str.contains(r"k|hill|constant|concentration")]
    assert check_correct_type(dataframe, float_cols, (float, int)) # both float and int are valid
    # Checks ratio
    ratio_cols = dataframe.columns[
        dataframe.columns.str.contains("ratio", case=False) &
        ~dataframe.columns.str.contains("concentration", case=False)] # anything that contains the substring ratio and not the substring concentration (ratio within concentration)
    for ratio_col in ratio_cols:
        dataframe[ratio_col] = dataframe[ratio_col].apply(define_ratio_from_string)
    # Checks regions
    regions_cols = dataframe.columns[dataframe.columns.str.contains("regions")]
    for regions_col in regions_cols:
        dataframe[regions_col] = dataframe[regions_col].apply(
            lambda x: define_region_list(x, num_regions))
    return dataframe

def extract_species_and_enzymes_from_reactions_data(
    enzymatic_reactions_df: pd.DataFrame, spontaneous_reactions_df: pd.DataFrame):
    """ Returns a numpy array with all the names of the involved species in the reactions defined 
    and a numpy array with all the names of the involved enzymes in the reactions defined."""
    unique_species = pd.concat([
        enzymatic_reactions_df['start_species'], enzymatic_reactions_df['end_species'],
        spontaneous_reactions_df['start_species'], spontaneous_reactions_df['end_species']
    ]).unique()
    unique_enzymes = enzymatic_reactions_df['enzyme'].unique()
    return unique_species, unique_enzymes

def extract_species_and_enzymes_from_characteristics_data(
        species_df: pd.DataFrame, enzymes_df: pd.DataFrame):
    """ Returns a numpy array of all the species defined in the species data and that 
    no species is defined more than once.
    """
    species = species_df['name']
    assert checks_lack_of_repetitions(species)

    enzymes = enzymes_df["name"]
    assert checks_lack_of_repetitions(enzymes)
    
    return list(species), list(enzymes)

def check_validity_reaction_network_info(case_directory: str, csv_file_names: list, num_regions: int):
    """
    case_directory: relative or absolute path to caseNNN
    csv_file_names: names of csv files without .csv

    Looks at the .csv files in the case_directory and checks the validity of the .csv files.
    True if all good, raises error if not.
    The .csv files are valid if:
    0. None of the .csv files has empty cells.
    1. Every species mentioned in the reaction files (enzymaticReactions.csv and spontaneousReactions.csv)
       is defined (only once) in species.csv. Inform of whether some species is defined that does not partake
       in any reaction.
    2. Every enzyme mentioned in the reaction files (enzymaticReactions.csv and spontaneousReactions.csv)
       is defined (only once) in enzymes.csv. Inform of whether some enzyme is defined that does not catalyze
       any reaction.
    If this is done, the dataframes are stored as pickled dataframes
    """
    # Create dataframes from .csv files (checks whether rows are duplicated)
    dataframes = {
        name: create_pandas_dataframe_from_csv_file(
        os.path.join(case_directory, f"{name}.csv"), num_regions)
        for name in csv_file_names
    }
    # Step 0
    for dataframe_name, dataframe_object in dataframes.items():
        assert no_empty_cells(dataframe_object), f"dataframe {dataframe_name} has empty cells"

    reaction_species, reaction_enzymes = extract_species_and_enzymes_from_reactions_data(dataframes["enzymatic_reactions"], dataframes["spontaneous_reactions"])
    species, enzymes = extract_species_and_enzymes_from_characteristics_data(dataframes["species"], dataframes["enzymes"])
    
    # Step 1
    missing_species_description = [
        individual_reaction_species for individual_reaction_species in reaction_species
        if individual_reaction_species not in species]
    if len(missing_species_description) == 0:
        print("All species that partake in reactions have their properties defined.")
    else:
        raise ImportError("The species", missing_species_description, "are not properly defined in species.csv")
    species_not_in_reaction = [individual_species for individual_species in species
        if individual_species not in reaction_species]
    if len(species_not_in_reaction) != 0:
        print("Species", species_not_in_reaction, "are defined but do not partake in any reaction.")
    # Step 2
    missing_enzymes_description = [
        individual_reaction_enzymes for individual_reaction_enzymes in reaction_enzymes
        if individual_reaction_enzymes not in enzymes]
    if len(missing_enzymes_description) == 0:
        print("All enzymes that catalyze reactions have their properties defined.")
    else:
        raise ImportError("The enzymes", missing_enzymes_description, "are not properly defined in ENZYMES.csv")
    enzymes_not_in_reaction = [individual_enzymes for individual_enzymes in enzymes
        if individual_enzymes not in reaction_enzymes]
    if len(enzymes_not_in_reaction) != 0:
        print("Enzymes", enzymes_not_in_reaction, "are defined but do not catalyze any reaction.")
    
    print(f"{case_directory} has good inputs. Pickling dataframes...")
    
    # Pickle dataframes separately
    for filename, dataframe in dataframes.items():
        pickle_dump_binary(
            os.path.join(case_directory, f".pickled_dataframe_{filename}"), dataframe)
    print(f"{case_directory} is good to go!")

if __name__ == "__main__":
    folder_to_check_validity = sys.argv[1]
    reaction_network_info_file_names = load_json("src/_template_reaction_network.json").keys()
    system_geometry_dict = load_json(
        os.path.join(folder_to_check_validity, ".expanded_system_geometry.json"))
    num_regions = system_geometry_dict["geometry_config"]["num_regions"]

    check_validity_reaction_network_info(folder_to_check_validity, reaction_network_info_file_names, num_regions)