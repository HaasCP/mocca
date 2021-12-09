from dataclasses import dataclass
from typing import Optional, Any, List
# TODO: Replace Any by DadData


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
    dataset : Any  # DADData parent of peak
    idx : int

    def __eq__(self, other):
        if not isinstance(other, BasePeak):
            # don't attempt to compare against unrelated types
            raise ValueError("Both peaks must be of the same type!")
        return (self.left == other.left and
                self.right == other.right and
                self.maximum == other.maximum and
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
class PreprocessedPeak(IntegratedPeak):
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
    dataset : Any  # DADData parent of peak
    idx : int
    saturation : bool
    pure : bool
    compound_id : Optional[str]
    integral : float
    concentration : Optional[float]

    def __eq__(self, other):
        if not isinstance(other, ProcessedPeak):
            # don't attempt to compare against unrelated types
            raise ValueError("Both peaks must be of the same type!")
        return (self.left == other.left and
                self.right == other.right and
                self.maximum == other.maximum and
                self.dataset == other.dataset)
