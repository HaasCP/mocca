from dataclasses import dataclass
from typing import Optional, List

import mocca.dad_data.models
import numpy as np

@dataclass(frozen=True)
class BasePeak():
    """
    Base peak class.
    """
    left : int
    right : int
    maximum : int


@dataclass(frozen=True)
class PickedPeak(BasePeak):
    """
    Class for picked peaks out of DAD data. Also valid for expanded peaks.
    """
    # https://www.python.org/dev/peps/pep-0484/#forward-references
    dataset : 'mocca.dad_data.models.CompoundData'  # DADData parent of peak
    idx : int

    def __eq__(self, other):
        if not isinstance(other, BasePeak):
            # don't attempt to compare against unrelated types
            raise ValueError("Both peaks must be of the same type!")
        return (self.maximum == other.maximum and
                self.dataset == other.dataset)


@dataclass(frozen=True, eq=False)
class CheckedPeak(PickedPeak):
    """
    Class for peaks checked with regard to saturation and purity.
    """
    saturation : bool
    pure : bool


@dataclass(frozen=True, eq=False)
class IntegratedPeak(CheckedPeak):
    """
    Class for integrated peaks.
    """
    integral : float


@dataclass(frozen=True, eq=False)
class CorrectedPeak(IntegratedPeak):
    """
    Class for peaks with added retention time offset. From this class on,
    retention times in the peaks are already corrected. This means, that accessing
    data from the dataset attribute require prior un-offsetting.
    """
    offset : int


@dataclass(frozen=True, eq=False)
class PreprocessedPeak(CorrectedPeak):
    """
    Class for preprocessed peaks containing a list of possible component matches
    in the attribute compound_id.
    """
    matches : List[dict]


@dataclass(frozen=True)
class ProcessedPeak():
    """
    Class of fully processed peaks ready to be put in the peak databse.
    """
    left : int
    right : int
    maximum : int
    dataset : 'mocca.dad_data.models.CompoundData'
    idx : int
    saturation : bool
    pure : bool
    integral : float
    offset : int
    compound_id : Optional[str] = None
    concentration : Optional[float] = None
    is_compound : bool = False

    def __eq__(self, other):
        if not isinstance(other, ProcessedPeak):
            # don't attempt to compare against unrelated types
            raise ValueError("Both peaks must be of the same type!")
        return (self.maximum + self.offset == other.maximum + other.offset and
                self.dataset == other.dataset)

@dataclass(frozen=True)
class ParafacPeak(BasePeak):
    integral : float
    spectra : np.ndarray
    elution : np.ndarray
    offset : int
