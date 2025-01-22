import gut_metabolites_labels as lb
import serum_metabolites_labels as sb
import shared_types as ty
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def subset_start_end_diagnosis(
    metabolites: str,
    diagnoses: str,
    gender: str,
    metabolite_group: ty.MetaboliteGroup | list[str],
    start_diagnosis: ty.Diagnosis,
    end_diagnosis: ty.Diagnosis,
) -> pd.DataFrame:
    """
    Returns a subset Dataframe based on starting and final diagnosis.
    """

    fatty_acid_data = _extract_metabolites_groups(
        metabolites_file=metabolites,
        diagnoses_file=diagnoses,
        gender_file=gender,
        metabolite_group=metabolite_group,
    )

    subset_end_diagnosis = _subset_baseline_from_enddx(
        data=fatty_acid_data, start_diagnosis=start_diagnosis
    )

    match end_diagnosis:
        case ty.Diagnosis.Normal:
            return subset_end_diagnosis.start_diagnosis.loc[
                subset_end_diagnosis.start_diagnosis["RID"].isin(
                    subset_end_diagnosis.start_to_normal
                )
            ]
        case ty.Diagnosis.MCI:
            return subset_end_diagnosis.start_diagnosis.loc[
                subset_end_diagnosis.start_diagnosis["RID"].isin(
                    subset_end_diagnosis.start_to_mci
                )
            ]
        case ty.Diagnosis.AD:
            return subset_end_diagnosis.start_diagnosis.loc[
                subset_end_diagnosis.start_diagnosis["RID"].isin(
                    subset_end_diagnosis.start_to_dementia
                )
            ]
        case ty.Diagnosis.All:
            return subset_end_diagnosis.start_diagnosis


def add_enddx_total_to_baseline(
    metabolites: str,
    diagnoses: str,
    gender: str,
    metabolite_group: ty.MetaboliteGroup | list[str],
    sample_source: ty.Sample,
    start_diagnosis: ty.Diagnosis,
    end_diagnosis: ty.Diagnosis = ty.Diagnosis.All,
) -> pd.DataFrame:
    """
    Adds the corresponding final diagnosis and total metabolites columns to baseline data.
    """
    assert start_diagnosis != ty.Diagnosis.All
    fatty_acid_data = _extract_metabolites_groups(
        metabolites_file=metabolites,
        diagnoses_file=diagnoses,
        gender_file=gender,
        metabolite_group=metabolite_group,
    )

    subset_end_diagnosis = _subset_baseline_from_enddx(
        data=fatty_acid_data, start_diagnosis=start_diagnosis
    )

    start_to_normal = subset_end_diagnosis.start_diagnosis.loc[
        subset_end_diagnosis.start_diagnosis["RID"].isin(
            subset_end_diagnosis.start_to_normal
        )
    ]
    baseline_start_normal = start_to_normal.loc[start_to_normal["VISCODE2"] == "bl"]
    baseline_start_normal.insert(loc=0, column="END_DIAGNOSIS", value="Normal")

    start_to_mci = subset_end_diagnosis.start_diagnosis.loc[
        subset_end_diagnosis.start_diagnosis["RID"].isin(
            subset_end_diagnosis.start_to_mci
        )
    ]
    baseline_start_mci = start_to_mci.loc[start_to_mci["VISCODE2"] == "bl"]
    baseline_start_mci.insert(loc=0, column="END_DIAGNOSIS", value="MCI")

    start_to_ad = subset_end_diagnosis.start_diagnosis.loc[
        subset_end_diagnosis.start_diagnosis["RID"].isin(
            subset_end_diagnosis.start_to_dementia
        )
    ]

    baseline_start_ad = start_to_ad.loc[start_to_ad["VISCODE2"] == "bl"]
    baseline_start_ad.insert(loc=0, column="END_DIAGNOSIS", value="AD")

    new_data = pd.concat(
        [baseline_start_normal, baseline_start_mci, baseline_start_ad], axis=0
    )

    variables = _name_variables(
        metabolite_group=metabolite_group, sample_source=sample_source
    )
    columns = list(variables.keys())[2:]

    new_data["Total_Metabolites"] = new_data.loc[
        :,
        columns,
    ].sum(axis=1)

    match end_diagnosis:
        case ty.Diagnosis.AD:
            new_data = new_data.loc[new_data["END_DIAGNOSIS"] == "AD"]
        case ty.Diagnosis.MCI:
            new_data = new_data.loc[new_data["END_DIAGNOSIS"] == "MCI"]
        case ty.Diagnosis.Normal:
            new_data = new_data.loc[new_data["END_DIAGNOSIS"] == "Normal"]

    return new_data


def display_pairwise_results_on_plot(
    ax, pairs: list, pvals: list, corrected_pvals, significance_level
):
    """Display pairwise comparison results on the plot."""
    pairwise_text_lines = []
    for (group1, group2), p, cp in zip(pairs, pvals, corrected_pvals):
        significance = "*" if cp < significance_level else ""
        pairwise_text_lines.append(
            f"  {group1} vs {group2}: adj p = {cp:.4f}{significance}"
        )

    bonferroni_text_lines = [f"Bonferroni Adj alpha = {significance_level:.4f}"]

    title_text = "\n".join(pairwise_text_lines + bonferroni_text_lines)

    ax.set_title(
        title_text,
        fontsize=10,
        verticalalignment="center",
        horizontalalignment="center",
    )

    ax.set_title(title_text, fontsize=10, ha="center")


def get_color_config(color_config: ty.FilterColor):
    """
    Returns the corresponding color configuration.
    """
    match color_config:
        case ty.FilterColor.Sex:
            return {"Male": "orchid", "Female": "darkorange"}

        case ty.FilterColor.Diagnosis:
            return {
                "MCI": "steelblue",
                "AD": "salmon",
                "Normal": "mediumseagreen",
            }


def plot_sunburst_dist(
    data: pd.DataFrame,
    path: list[str],
    color_sequence: list[str] = ["steelblue", "salmon", "mediumseagrean"],
):
    """
    Returns a sunburst plot based on the declared path and data.
    """
    fig = go.Figure(
        px.sunburst(data, path=path, color_discrete_sequence=color_sequence)
    )

    fig.update_layout(width=800, height=800)
    fig.update_traces(textfont_size=20)
    fig.update_traces(textinfo="label+value+percent parent")
    fig.update_layout(margin=dict(t=10, l=10, r=10, b=10))
    fig.show()


def _extract_and_clean_data(
    data: pd.DataFrame, columns: list[str], right_data: pd.DataFrame
) -> pd.DataFrame:
    """
    Returns a full dataframe with demographic data, duplicates imputed with mean.
    """

    return (
        (((data[columns]).groupby(["RID", "VISCODE2"])).mean())
        # impute duplicates with mean
        .merge(
            right=right_data,
            how="left",
            left_on=["RID", "VISCODE2"],
            right_on=["RID", "VISCODE2"],
        )
    )


def _get_metabolite_group(metabolite_group: ty.MetaboliteGroup):
    """
    Returns a dictionary of column names for metabolite groups.
    """
    match metabolite_group:
        case ty.MetaboliteGroup.ShortFattyAcids:
            return lb.short_chained
        case ty.MetaboliteGroup.MediumFattyAcids:
            return lb.medium_chained
        case ty.MetaboliteGroup.LongFattyAcids:
            return lb.long_chained
        case ty.MetaboliteGroup.AminoAcids:
            return lb.amino_acids
        case ty.MetaboliteGroup.OrganicAcids:
            return lb.organic_acids
        case ty.MetaboliteGroup.IndoleCompounds:
            return lb.indole_compounds
        case ty.MetaboliteGroup.BileAcids:
            return lb.bile_acids
        case ty.MetaboliteGroup.OxidizedBileAcids:
            return lb.oxidized_bile_acids
        case ty.MetaboliteGroup.GlycineConjugatedBileAcids:
            return lb.glycine_conjugated
        case ty.MetaboliteGroup.TaurineConjugatedBileAcids:
            return lb.taurine_conjugated
        case ty.MetaboliteGroup.IsomericRareBileAcids:
            return lb.isomeric_and_rare
        case ty.MetaboliteGroup.Others:
            return lb.others
        case ty.MetaboliteGroup.All:
            return lb.all_metabolites
        case ty.MetaboliteGroup.SerumCholesterol:
            return sb.cholesterol
        case ty.MetaboliteGroup.SerumTriglycerides:
            return sb.triglycerides
        case ty.MetaboliteGroup.SerumPhospholipids:
            return sb.phospholipids
        case ty.MetaboliteGroup.SerumCholesterylEsters:
            return sb.cholesteryl_esters
        case ty.MetaboliteGroup.SerumFreeCholesterol:
            return sb.free_cholesterol
        case ty.MetaboliteGroup.SerumLipids:
            return sb.lipids
        case ty.MetaboliteGroup.SerumLipoproteins:
            return sb.lipoproteins
        case ty.MetaboliteGroup.SerumOtherLipids:
            return sb.other_lipids
        case ty.MetaboliteGroup.SerumAlipoproteins:
            return sb.alipoproteins
        case ty.MetaboliteGroup.SerumFattyAcids:
            return sb.fatty_acids
        case ty.MetaboliteGroup.SerumAminoAcids:
            return sb.amino_acids
        case ty.MetaboliteGroup.SerumGlycolysis:
            return sb.glycolysis_related
        case ty.MetaboliteGroup.SerumKetones:
            return sb.ketone_bodies
        case ty.MetaboliteGroup.SerumBalanceInflamm:
            return sb.fluid_balance_inflammation
        case ty.MetaboliteGroup.SerumXXLVLDL:
            return sb.extremely_large_VLDL
        case ty.MetaboliteGroup.SerumXLVLDL:
            return sb.very_large_VLDL
        case ty.MetaboliteGroup.SerumLVLDL:
            return sb.large_VLDL
        case ty.MetaboliteGroup.SerumMVLDL:
            return sb.medium_VLDL
        case ty.MetaboliteGroup.SerumSVLDL:
            return sb.small_VLDL
        case ty.MetaboliteGroup.SerumXSVLDL:
            return sb.very_small_VLDL
        case ty.MetaboliteGroup.SerumIDL:
            return sb.idl
        case ty.MetaboliteGroup.SerumLLDL:
            return sb.large_LDL
        case ty.MetaboliteGroup.SerumMLDL:
            return sb.medium_LDL
        case ty.MetaboliteGroup.SerumSLDL:
            return sb.small_LDL
        case ty.MetaboliteGroup.SerumXLHDL:
            return sb.very_large_HDL
        case ty.MetaboliteGroup.SerumLHDL:
            return sb.large_HDL
        case ty.MetaboliteGroup.SerumMHDL:
            return sb.medium_HDL
        case ty.MetaboliteGroup.SerumSHDL:
            return sb.small_HDL
        case ty.MetaboliteGroup.SerumAll:
            return sb.all_metabolites


def _extract_metabolites_groups(
    metabolites_file: str,
    diagnoses_file: str,
    gender_file: str,
    metabolite_group: ty.MetaboliteGroup | list[str],
) -> pd.DataFrame:
    """Returns a subset of the metabolites dataframe based on fatty acid chain classification\n"""
    """The duplicate entries were replaced by the mean of the results"""

    metabolites = pd.read_csv(
        metabolites_file,
        low_memory=False,
    )
    diagnoses = (
        pd.read_csv(
            diagnoses_file,
            low_memory=False,
        )
    )[["RID", "VISCODE2", "DIAGNOSIS"]]
    gender_and_age = (
        pd.read_csv(
            gender_file,
            low_memory=False,
        )
    )[["RID", "VISCODE", "PTGENDER", "AGE"]].rename(columns={"VISCODE": "VISCODE2"})
    dataframe = metabolites.copy()
    merged_gender_age_and_dx = diagnoses.merge(
        right=gender_and_age, how="left", on=["RID", "VISCODE2"]
    )

    if isinstance(metabolite_group, list):
        return _extract_and_clean_data(
            data=dataframe,
            right_data=merged_gender_age_and_dx,
            columns=(["RID", "VISCODE2"] + metabolite_group),
        )

    else:
        return _extract_and_clean_data(
            data=dataframe,
            right_data=merged_gender_age_and_dx,
            columns=list(
                (_get_metabolite_group(metabolite_group=metabolite_group)).keys()
            ),
        )


def _subset_baseline_from_enddx(data: pd.DataFrame, start_diagnosis: ty.Diagnosis):
    baseline_data = data.loc[data["VISCODE2"] == "bl"]
    start = data.loc[
        data["RID"].isin(
            baseline_data.loc[baseline_data["DIAGNOSIS"] == start_diagnosis.value][
                "RID"
            ]
        )
    ]
    normie = set((start.loc[start["DIAGNOSIS"] == 1.0])["RID"])
    mci = set((start.loc[start["DIAGNOSIS"] == 2.0])["RID"])
    ad = set((start.loc[start["DIAGNOSIS"] == 3.0])["RID"])
    start_to_dementia = ad
    start_to_mci = mci - ad
    start_to_normal = normie - mci - ad
    return ty.SubsetEndDiagnosis(
        start, start_to_normal, start_to_mci, start_to_dementia
    )


def _name_variables(
    metabolite_group: ty.MetaboliteGroup | list[str], sample_source: ty.Sample
):
    if isinstance(metabolite_group, list):
        match sample_source:
            case ty.Sample.Gut:
                assert set(metabolite_group) <= set(lb.all_metabolites.keys())
                variable_names = ["RID", "VISCODE2"] + metabolite_group
                return {
                    metabolite: lb.all_metabolites[metabolite]
                    for metabolite in variable_names
                }
            case ty.Sample.Serum:
                assert set(metabolite_group) <= set(sb.all_metabolites.keys())
                variable_names = ["RID", "VISCODE2"] + metabolite_group
                return {
                    metabolite: sb.all_metabolites[metabolite]
                    for metabolite in variable_names
                }

    else:
        return _get_metabolite_group(metabolite_group=metabolite_group)
