#!/usr/bin/env python

import argparse
import csv
import ete3

ranks = ["-", "D", "K", "P", "C", "O", "F", "G", "S"]

if __name__=='__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', help="input combined bracken report", required=True)
    parser.add_argument('-o', help="output newick tree", required=True)

    args = parser.parse_args()

    nodes = {}

    tree = ete3.Tree()
    root = tree.add_child(name='root')

    nodes[0] = root

    with open(args.i) as handle:
        for line in csv.reader(handle, delimiter='\t'):
            fraction, cumulative, count, order, tax_id, taxa_entry = line

            if not order in ranks:
                continue

            if order == 'D': #new domain, start with new branch
                nodes = {0: root}

            taxa_name = taxa_entry.strip()
            index = ranks.index(order)
            #Certain levels don't exist everywhere, e.g. no K in Bacteria
            parent = nodes[index-1] if index-1 in nodes.keys() else nodes[index-2]

            child = parent.add_child(name=tax_id)
            nodes[index] = child
    
    tree.write(outfile=args.o, format=3)
