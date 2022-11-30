def xr_to_csv(dset, csv_stream):
    df = dset.to_dataframe()
    df.to_csv(csv_stream, line_terminator='\n')


def csv_to_xr(csv_stream, index_cols):
    import pandas as pd
    import xarray as xr

    df = pd.read_csv(csv_stream, index_col=index_cols)
    dset = xr.Dataset.from_dataframe(df)
    return dset
