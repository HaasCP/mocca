import pandas as pd


def read_csv(csv_path, encoding):
    """
    Reads the csv files with given encoding.  
    Parameters
    ----------
    csv_path : str
        The path to the csv file.

    encoding: str
        The encoding of the experimental data.

    Returns
    -------
    df : pandas.DataFrame
        First column is time, the following columns obtain the absorbance 
        values at the given detection wavelength in the column name.
    """

    with open(csv_path, 'r', encoding=encoding) as f:
        df = pd.read_csv(f)
    return df
