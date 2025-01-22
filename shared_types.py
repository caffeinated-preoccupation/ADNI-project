from dataclasses import dataclass
from enum import Enum
import pandas as pd


class BioSex(Enum):
    Male = 1.0
    Female = 2.0
    All = 3.0


class Diagnosis(Enum):
    Normal = 1.0
    MCI = 2.0
    AD = 3.0
    All = 4.0


class Sample(Enum):
    Gut = 1.0
    Serum = 2.0


class MetaboliteGroup(Enum):
    ShortFattyAcids = 1.0
    MediumFattyAcids = 2.0
    LongFattyAcids = 3.0
    AminoAcids = 4.0
    OrganicAcids = 5.0
    IndoleCompounds = 6.0
    BileAcids = 7.0
    OxidizedBileAcids = 8.0
    GlycineConjugatedBileAcids = 9.0
    TaurineConjugatedBileAcids = 10.0
    IsomericRareBileAcids = 11.0
    Others = 12.0
    All = 13.0
    SerumCholesterol = 14.0
    SerumTriglycerides = 15.0
    SerumPhospholipids = 16.0
    SerumCholesterylEsters = 17.0
    SerumFreeCholesterol = 18.0
    SerumLipids = 19.0
    SerumLipoproteins = 20.0
    SerumOtherLipids = 21.0
    SerumAlipoproteins = 22.0
    SerumFattyAcids = 23.0
    SerumAminoAcids = 24.0
    SerumGlycolysis = 25.0
    SerumKetones = 26.0
    SerumBalanceInflamm = 27.0
    SerumXXLVLDL = 28.0
    SerumXLVLDL = 29.0
    SerumLVLDL = 30.0
    SerumMVLDL = 31.0
    SerumSVLDL = 32.0
    SerumXSVLDL = 33.0
    SerumIDL = 34.0
    SerumLLDL = 35.0
    SerumMLDL = 36.0
    SerumSLDL = 37.0
    SerumXLHDL = 38.0
    SerumLHDL = 39.0
    SerumMHDL = 40.0
    SerumSHDL = 41.0
    SerumAll = 42.0


class FilterColor(Enum):
    Sex = 1.0
    Diagnosis = 2.0


@dataclass(frozen=True, slots=True)
class SubsetEndDiagnosis:
    start_diagnosis: pd.DataFrame
    start_to_normal: set
    start_to_mci: set
    start_to_dementia: set
