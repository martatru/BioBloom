#!/usr/bin/env python3
import os, sys, glob

def main():
    folder = input("Enter the path to the folder with CSV files: ").strip().strip('"')
    if not folder or not os.path.isdir(folder):
        print(f"Not a directory: {folder}")
        sys.exit(1)

    output_path = os.path.join(folder, "merged.csv")
    csv_files = sorted(glob.glob(os.path.join(folder, "*.csv")))
    csv_files = [p for p in csv_files if os.path.abspath(p) != os.path.abspath(output_path)]
    if not csv_files:
        print("No .csv files found in the specified folder.")
        sys.exit(0)

    wrote_header = False
    lines_written = 0
    files_merged = 0

    with open(output_path, "w", encoding="utf-8", newline="") as out_f:
        for path in csv_files:
            with open(path, "r", encoding="utf-8-sig", newline="") as f:
                first_line = f.readline()
                if first_line == "":
                    continue
                if not wrote_header:
                    out_f.write(first_line)
                    lines_written += 1
                    for line in f:
                        out_f.write(line)
                        lines_written += 1
                    wrote_header = True
                else:
                    for line in f:
                        out_f.write(line)
                        lines_written += 1
                files_merged += 1

    if not wrote_header:
        print("All CSV files were empty. Created an empty merged.csv.")
    else:
        print(f"Done. Merged {files_merged} file(s) into: {output_path}")
        print(f"Total number of lines in output file (including header): {lines_written}")

if __name__ == "__main__":
    main()
