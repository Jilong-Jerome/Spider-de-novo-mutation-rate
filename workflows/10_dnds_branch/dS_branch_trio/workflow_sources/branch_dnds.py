#!/usr/bin/env python3
"""
branch_dnds.py - Extract per-branch dN/dS from paml2tab.py output

For trio tree ((SAR #1,PAC #2) #4,BIC #3), paml2tab.py produces rows for:
  SAR, PAC, BIC, SAR_PAC (internal node)

This script reads the tab file and outputs per-branch dN, dS, dN/dS.

Usage:
    python3 branch_dnds.py input.tab output_branch.tab
"""

import sys


def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    data = {}
    with open(input_file) as f:
        header = next(f)
        for line in f:
            fields = line.strip().split('\t')
            nodename = fields[0]
            dN = float(fields[1])
            dS = float(fields[2])
            rep_id = fields[3]
            data[nodename] = {'dN': dN, 'dS': dS, 'id': rep_id}

    with open(output_file, 'w') as out:
        out.write("id\tbranch\tdN\tdS\tdNdS\n")
        for nodename, vals in data.items():
            dN = vals['dN']
            dS = vals['dS']
            rep_id = vals['id']
            if dS > 0:
                dNdS = dN / dS
            else:
                dNdS = 'NA'
            out.write(f"{rep_id}\t{nodename}\t{dN}\t{dS}\t{dNdS}\n")


if __name__ == '__main__':
    main()
