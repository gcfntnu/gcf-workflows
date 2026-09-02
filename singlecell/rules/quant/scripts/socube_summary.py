import sys

import pandas as pd

if __name__ == "__main__":
    df = pd.read_csv(sys.argv[1], index_col=0)
    df.index.name = "Barcode"
    df.rename(columns={"predict_score": "doublet_score", "predict_type": "doublet"}, inplace=True)
    df = df.iloc[:,[1,0]]
    df.reset_index().to_csv(sys.stdout, sep="\t", index=False)
