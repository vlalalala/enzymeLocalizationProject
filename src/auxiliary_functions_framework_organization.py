import ast

def get_dict_with_correct_key_types_from_json_file(
        imported_dict_from_json, species_lookup_dict) -> dict:
    """ Takes a nested dict imported from a json file with the structure
    [string key][string key][species.name string key]
    and returns the same nested dict but with the last key
    as [int key][int key][species object] using the dictionary species_lookup_dict, which maps
    the species name to the species object
    (also makes the integer keys that got converted to strings when originally saving the
    dict into a json file into integers again).

    Used for saved concentrations and point_ids.
    """
    dict_in_correct_concentrations_format = {
        int(region_idx): {
            int(mesh_point_idx): {
                species_lookup_dict[species_name]: data
                for species_name, data in mesh_point_info.items()}
            for mesh_point_idx, mesh_point_info in region_info.items()}
        for region_idx, region_info in imported_dict_from_json.items()
    }
    return dict_in_correct_concentrations_format

def get_correct_reverse_point_ids_dict(imported_reverse_point_ids_dict,
        species_lookup_dict
    ):
    dict_inc_correct_format = {
        int(k): [*v[:-1], species_lookup_dict[v[-1]]]
        for k, v in imported_reverse_point_ids_dict.items()
    }
    return dict_inc_correct_format

def get_correct_neighbors_dict(imported_neighbors_dict):
    """From the json file, the neighbors tuple gets saved as a list. Convert back.
    """
    new_dict = {}
    for k, v in imported_neighbors_dict.items():
        # Parse the string "[0, 0]" into a real list [0, 0]
        key_as_list = ast.literal_eval(k)
        # Convert it to a tuple
        key_as_tuple = tuple(key_as_list)
        new_dict[key_as_tuple] = v
    return new_dict

