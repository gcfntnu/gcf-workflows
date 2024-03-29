import sys

import pandas as pd

if __name__ == "__main__":
    df = pd.read_csv(sys.argv[1], index_col=0)
    df.index = ["{}-1".format(i.split("-")[0]) for i in df.index]
    df.index.name = "Barcode"
    df.reset_index().to_csv(sys.stdout, sep="\t", index=False)
