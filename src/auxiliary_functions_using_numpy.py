import sys
import numpy as np
import csv
import os


def sparse_txt_to_csv(input_txt_path, output_csv_path):
    """
    Reads a sparse matrix stored as:
        row col value
    and writes a dense CSV file where missing elements are 0.
    """

    rows = []
    cols = []
    values = []

    # Read the sparse file
    with open(input_txt_path, 'r') as f:
        next(f)  # skip header if present (row col value)
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            r, c, v = int(parts[0]), int(parts[1]), float(parts[2])
            rows.append(r)
            cols.append(c)
            values.append(v)

    # Determine matrix size
    n_rows = max(rows) + 1
    n_cols = max(cols) + 1

    # Create dense matrix initialized with zeros
    matrix = np.zeros((n_rows, n_cols))

    # Fill known values
    for r, c, v in zip(rows, cols, values):
        matrix[r, c] = v

    # Save to CSV
    with open(output_csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(matrix)

    print(f"Dense matrix saved to {output_csv_path}")



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sparse_to_csv.py input.txt output.csv")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    sparse_txt_to_csv(input_file, output_file)