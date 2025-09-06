import polars as pl

metadata = pl.read_csv("data/raw/Metadata.csv")

print(metadata.head())