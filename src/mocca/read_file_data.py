# flake8: noqa
import pandas as pd
import numpy as np


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

def absorbance_to_array(df):
    """
    Generates a 2D absorbance array of the absorbance values.

    Parameters
    ----------
    df : pandas.DataFrame
        Columns:
            time: Chromatogram time
            wavelength: Detection wavelength
            absorbance: absorbance value.

    Returns
    -------
    absorbance_array : numpy.ndarray
        Absorbance values with number of wavelengths in the first 
        and number of recorded times in the second dimension.

    """
    absorbance_array = df.absorbance.to_numpy().\
        reshape(df.wavelength.nunique(), df.time.nunique())
    return absorbance_array

def apply_filter(dataframe):
    """
    Filters absorbance data of tidy 3D DAD dataframes to remove noise
    and background systematic error.
    Parameters
    ----------
    dataframe : pandas.DataFrame
        Columns:
            time: Chromatogram time
            wavelength: Detection wavelength
            absorbance: absorbance value.

    Returns
    -------
    df : pandas.DataFrame
        Columns:
            time: Chromatogram time
            wavelength: Detection wavelength
            absorbance: absorbance value.
    """

    df = dataframe.copy()
    df['absorbance'] = df.groupby('time')['absorbance']\
        .rolling(window=5, center=True).mean().reset_index(0,drop=True)
    df = df.dropna().reset_index(0,drop=True)
    n_times = len(df.time.unique())
    wls = df.wavelength.unique()
    df['absorbance'] = df.absorbance - df[df['wavelength'] == wls.max()]\
        .absorbance.iloc[np.tile(np.arange(n_times), len(wls))]\
        .reset_index(0,drop=True)
    return df
