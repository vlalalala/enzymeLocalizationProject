import sys
import os
from itertools import product
from pathlib import Path
import pandas as pd
import numpy as np
from auxiliary_functions_using_standard_library import is_int_value
import numpy as np
from auxiliary_functions import read_yaml_file, dump_in_yaml_file
import pandas as pd
from itertools import product
from pathlib import Path
import yaml
import itertools

def parse_value_from_yaml(value):
    """Takes a value from a yaml and returns Python version of the value.
    The value can either be a dict (the value is a nested dict) or 
    a list.
    """
    if isinstance(value, dict): # Go one level further
        return {k: parse_value_from_yaml(val) for k, val in value.items()}
    if not isinstance(value, list):
        raise ValueError("Something went wrong parsing the yaml file. No value should not be within a list.")
    return parse_value(value)

def get_yaml_contents_as_dict(file_path):
    """Returns contents of yaml file as a python dictionary
    with necessary linspace & range conversions.
    """
    raw = read_yaml_file(file_path)
    expanded = {k: parse_value_from_yaml(v) for k, v in raw.items()}
    return expanded

def parse_value(raw):
    """ Used to parse a raw value into a python list of options. The value must be
    either a list (value of final nested dictionary) or a (child) dictionary.
    if value is ["range:0:3:1", 2], return is [0,1,2,3,2]
    if value is ["range:0:3:1", [4,5]], return is [0,1,2,3,[4,5]] (undesirable)
    if value is ["[range:0:3:1]", [4,5]], return is [[0,1,2,3], [4,5]] (desirable)
    Works analogously with linspace
    e.g. value is ["[linspace:1:4:3]", "[range:1:3:2]", [6,8]] -> [[1.0,2.5,4.0], [1,3], [6,8]]
    """
    if isinstance(raw, list): # case that raw comes from a yaml file
        items = raw
    elif isinstance(raw, str):
        s = raw.strip()
        loaded = yaml.safe_load(s)  # parses "[a, b]" into ["a","b"], "3" into 3, etc.
        items = loaded if isinstance(loaded, list) else [loaded]
    else:
        items = [raw]

    complete_values_list = []
    for item in items:
        if isinstance(item, str):
            string_encodes_list = False
            if item.startswith("[") and item.endswith("]"):
                string_encodes_list = True
                item = item[1:-1] # rename element to not have []
            if item.startswith("range:"):
                _, start, stop, step = item.split(":")
                if not (is_int_value(start) and is_int_value(stop) and is_int_value(step)):
                    raise ValueError(f"The range given through {item} is not valid.")
                # create a list with the range values
                final_elements_from_item = list(range(int(start), int(stop)+1, int(step)))
                if string_encodes_list:
                    complete_values_list.append(final_elements_from_item)
                    continue
                # if a list is not encoded,
                # append the individual elements from the range values to the complete list
                for final_element in final_elements_from_item:
                    complete_values_list.append(final_element)
            elif item.startswith("linspace"):
                is_list_mode = item.startswith("linspace_list:")
                if is_list_mode:
                    _, start, stop, num = item.split(":")
                else:
                    _, start, stop, num = item.split(":")
                if not is_int_value(num):
                    raise ValueError(f"The number of array elements given through {item} is not valid.")
                values = [float(x) for x in np.linspace(float(start), float(stop), int(num))]
                if is_list_mode:
                    values = [[v] for v in values]
                if string_encodes_list:
                    complete_values_list.append(values)
                else:
                    complete_values_list.extend(values)
            elif item.startswith("logspace:"):
                _, start, stop, num = item.split(":")
                if not is_int_value(num):
                    raise ValueError(f"The number of array elements given through {item} is not valid.")
                final_elements_from_item = list(
                    float(x) for x in np.logspace(float(start), float(stop), int(num)))
                if string_encodes_list:
                    complete_values_list.append(final_elements_from_item)
                    continue
                for x in final_elements_from_item:
                    complete_values_list.append(x)
            elif item.startswith("dictlinspace:"):
                """k1 is key 1, k2 is key 2, start is the initial value corresponding to k1,
                stop is the final value corresponding to k2, num is the number of elements in the list.
                Additional keys (k3, k4, ...) are assigned value 0.

                dictlinspace:0:1:2:3:0.0:1.0:5 → keys 0 and 1 get x/1-x, keys 2 and 3 get 0.0
                """
                parts = item.split(":")
                # parts[0] = "dictlinspace", parts[1] = k1, parts[2] = k2,
                # parts[-3] = start, parts[-2] = stop, parts[-1] = num
                # parts[3:-3] = optional extra keys with value 0
                if len(parts) < 6:
                    raise ValueError(f"Invalid dictlinspace item: {item}")
                k1, k2 = parts[1], parts[2]
                extra_keys = parts[3:-3]
                start, stop, num = parts[-3], parts[-2], parts[-1]
                if not is_int_value(num):
                    raise ValueError(f"Invalid num in {item}")
                values = np.linspace(float(start), float(stop), int(num))
                dicts = []
                for v in values:
                    d = {int(k1): float(v), int(k2): float(1 - v)}
                    for ek in extra_keys:
                        d[int(ek)] = 0.0
                    dicts.append(d)
                if string_encodes_list:
                    complete_values_list.append(dicts)
                else:
                    complete_values_list.extend(dicts)
            else:
                complete_values_list.append(item)
        else:
            complete_values_list.append(item)
    return complete_values_list

def scalar_to_csv_str(x):
    """
    Convert a selected scalar value to a CSV-safe string.
    """
    if x is None:
        return ""
    return str(x)

def write_one_csv_per_combination(
    input_csv: str,
    output_dir: str,
    max_files=None,
    single_option_columns=None,   # e.g. ["id", "name", "fixed_param"]
):  
    single_option_columns = set(single_option_columns or [])

    df = pd.read_csv(input_csv, dtype=str)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)


    # Collect expandable cells + validate restricted columns
    choice_cells = []
    for r in range(len(df)):
        for c in df.columns:
            opts = parse_value(df.at[r, c])

            # Enforce single-option restriction on specified columns
            if c in single_option_columns and len(opts) > 1:
                raise ValueError(
                    f"Column '{c}' must have exactly one option per cell, "
                    f"but row {r} has {len(opts)} options: {opts}"
                )

            if len(opts) > 1:
                choice_cells.append((r, c, opts))

    # If no choices, still write exactly one CSV
    if not choice_cells:
        out = df.copy()
        for r in range(len(out)):
            for c in out.columns:
                out.at[r, c] = scalar_to_csv_str(parse_value(out.at[r, c])[0])
        out.to_csv(out_dir / f"combo_000001.csv", index=False)
        return 1

    option_lists = [opts for (_, _, opts) in choice_cells]

    count = 0
    for combo in product(*option_lists):
        count += 1
        if max_files is not None and count > max_files:
            break

        out_df = df.copy(deep=True)

        # Fill selected values for choice-cells
        for (r, c, _opts), selected in zip(choice_cells, combo):
            out_df.at[r, c] = scalar_to_csv_str(selected)

        # Normalize all other cells to their single selected value
        for r in range(len(out_df)):
            for c in out_df.columns:
                if not any(r == rr and c == cc for rr, cc, _ in choice_cells):
                    out_df.at[r, c] = scalar_to_csv_str(parse_value(out_df.at[r, c])[0])

        out_df.to_csv(out_dir / f"combo_{count:06d}.csv", index=False)

    return count

def expand_nested_combinations(dictionary):
    """
    Take a nested dict where all leaves are lists,
    return a list of dicts where leaves are single values.
    """
    # 1. Recursively collect all leaf paths and their lists
    def collect_paths(prefix, obj, out):
        if isinstance(obj, dict):
            for k, v in obj.items():
                collect_paths(prefix + (k,), v, out)
        else:
            # leaf → must be a list of values
            out.append((prefix, obj))
    leaf_paths = []
    collect_paths((), dictionary, leaf_paths)
    # Verify leaves are lists
    for path, vals in leaf_paths:
        if not isinstance(vals, list):
            raise ValueError(f"Leaf at {path} is not a list: {vals}")
    # 2. Cartesian product of all leaf value lists
    all_value_lists = [vals for _, vals in leaf_paths]
    all_combos = list(itertools.product(*all_value_lists))
    # 3. Rebuild nested dict for each combination
    results = []
    for combo in all_combos:
        newdict = {}
        for (path, _), value in zip(leaf_paths, combo):
            # walk down path creating structure
            cur = newdict
            for key in path[:-1]:
                cur = cur.setdefault(key, {})
            cur[path[-1]] = value
        results.append(newdict)
    return results

def write_one_yaml_file_per_combination(
        input_yaml_file: str, output_dir: str):
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Create dictionaries with only one value
    dictionary = get_yaml_contents_as_dict(input_yaml_file)
    single_values_dict_list = expand_nested_combinations(dictionary)
    for single_values_dict_index, single_values_dict in enumerate(single_values_dict_list):
        dump_in_yaml_file(out_dir / f"combo_{single_values_dict_index:06d}.yaml", single_values_dict)

if __name__ == "__main__":
    folder_with_parameter_ranges = sys.argv[1]
    for filename in os.listdir(folder_with_parameter_ranges):
        name, extension = os.path.splitext(filename)
        if extension == '.csv':
            write_one_csv_per_combination(
                os.path.join(folder_with_parameter_ranges, filename),
                os.path.join(folder_with_parameter_ranges, f"options_{name}"),
                single_option_columns=["name", "enzyme", "start_species", "end_species"]
            )
        elif extension == '.yaml':
            write_one_yaml_file_per_combination(
                os.path.join(folder_with_parameter_ranges, filename),
                os.path.join(folder_with_parameter_ranges, f"options_{name}"),
            )

        

            

