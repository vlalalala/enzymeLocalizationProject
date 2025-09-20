import os
import re

# Create new case

# Step 1: List all subfolders in 'data/' that match 'caseNNN'
existing_cases = [
    int(re.search(r"case(\d+)", d).group(1))
    for d in os.listdir(".")
    if os.path.isdir(d) and re.match(r"case\d+$", d)
]
# Step 2: Find the next case number
if existing_cases:
    next_case_number = max(existing_cases) + 1
else:
    next_case_number = 0

# Step 3: Format the new folder name
next_case = f"case{next_case_number:03d}"
next_case_dir = f"data/{next_case}"

rule create_next_case:
    output:
        f"{next_case_dir}/enzymaticReactions.csv {next_case_dir}/spontaneousReactions.csv"
    shell:
        """
        mkdir -p {os.path.dirname(output)}
        echo "Generated in {output}" > {output}
        """