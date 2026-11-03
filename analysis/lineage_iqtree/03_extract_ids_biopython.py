#!/usr/bin/env python3
#                 ############################
#               ##                            ##
#             ##  ####                    ####  ##
#           ##   ######                  ######   ##
#         ##    ########                ########    ##
#       ##      ##########            ##########      ##
#     ##        ####  ####################  ####        ##
#   ##          #####  ##################  #####          ##
# ##            ################################            ##
#   ##            ############################            ##
#     ##              ####################              ##
#       ##             ##  ##########  ##             ##
#         ##           ##   ########   ##           ##
#           ##          #    ######    #          ##
#             ##              ####              ##
#               ##                            ##
#                 ############################
#
# ==========================================================
#  PROJECT : URGG catalog
#  SCRIPT  : 03_extract_ids_biopython.py  -  Extract selected sequences from an alignment
#  AUTHOR  : QRQ
#  INPUT   : MSA FASTA, ID list (arguments)
#  OUTPUT  : subset FASTA
# ==========================================================

import gzip
import sys

from Bio import SeqIO

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def main():
    src, id_list_path, out_path = sys.argv[1:4]

    with open(id_list_path) as f:
        wanted = {line.strip() for line in f if line.strip()}

    found = set()
    with open_maybe_gz(src) as fh, open(out_path, "w") as out:
        for rec in SeqIO.parse(fh, "fasta"):
            # tolerate RS_/GB_ prefixes GTDB-Tk sometimes adds to headers
            bare_id = rec.id.replace("RS_", "").replace("GB_", "")
            if rec.id in wanted or bare_id in wanted:
                SeqIO.write(rec, out, "fasta")
                found.add(rec.id if rec.id in wanted else bare_id)

    missing = wanted - found
    print(f"Matched {len(found)}/{len(wanted)} IDs.")
    if missing:
        print(f"Missing {len(missing)} IDs (first 10 shown): {sorted(missing)[:10]}")

if __name__ == "__main__":
    main()
