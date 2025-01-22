import math
from itertools import combinations

import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from scipy.stats import f_oneway
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests
import numpy as np
import shared_functions as fu
import shared_types as ty
import pandas as pd


def plot_hist_kde(
    metabolite_group: ty.MetaboliteGroup | list[str],
    hue: ty.FilterColor,
    start_diagnosis: ty.Diagnosis,
    sample_source: ty.Sample,
    end_diagnosis: ty.Diagnosis = ty.Diagnosis.All,
    biological_sex: ty.BioSex = ty.BioSex.All,
    color_config: ty.FilterColor = ty.FilterColor.Sex,
):
    """
    Plots the histogram & density plots of the variables.

    Arguments:
    - metabolite_group: This takes on a class or a list of column names from the original dataset.
    - hue: The column to group data by.
    - start_diagnosis: The starting diagnosis filter.
    - end_diagnosis: The ending diagnosis filter (default is Diagnosis.All).
    - biological_sex: The gender filter (default is BioSex.All).
    - legend_title: Custom title for the legend.
    """
    assert start_diagnosis != ty.Diagnosis.All

    data, hue_value, variables, hist_color = _prepare_data_and_variables(
        metabolite_group=metabolite_group,
        hue=hue,
        start_diagnosis=start_diagnosis,
        sample_source=sample_source,
        end_diagnosis=end_diagnosis,
        biological_sex=biological_sex,
        color_config=color_config,
    )

    columns = 4
    rows = math.ceil(len(variables) / columns)  # Determine the required number of rows
    fig, axes = plt.subplots(nrows=rows, ncols=columns, figsize=(6 * columns, 6 * rows))
    axes = axes.ravel()

    for i, (variable, title) in enumerate(list(variables.items())[2:]):
        model_formula = f"{variable} ~ C({hue_value})"
        model = ols(model_formula, data=data).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)  # Perform ANOVA
        p_value = round(anova_table.loc[f"C({hue_value})", "PR(>F)"], 4)

        sns.histplot(
            data=data,
            x=variable,
            hue=hue_value,
            ax=axes[i],
            kde=True,
            stat="density",
            common_norm=False,
            element="step",
            palette=hist_color,
        )
        axes[i].set_title(f"ANOVA p-value = {p_value:.4f}")
        axes[i].set_ylabel("Density")
        axes[i].set_xlabel(title)

        if p_value < 0.05 and hue == ty.FilterColor.Diagnosis:
            group_names = data[hue_value].unique()
            pairwise_pvalues = []
            pairs = list(combinations(group_names, 2))

            for group1, group2 in pairs:
                group1_data = data[data[hue_value] == group1][variable]
                group2_data = data[data[hue_value] == group2][variable]
                _, p = f_oneway(group1_data, group2_data)
                pairwise_pvalues.append(p)

            # Bonferroni correction
            _, corrected_pvalues, _, _ = multipletests(
                pairwise_pvalues, method="bonferroni"
            )
            significance_level = 0.05 / len(pairwise_pvalues)

            # Update footnote with adjusted significance level
            fu.display_pairwise_results_on_plot(
                axes[i], pairs, pairwise_pvalues, corrected_pvalues, significance_level
            )

    # Hide any unused subplots
    for j in range(i + 1, rows * columns):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.show()


def display_effect_sizes_and_cohens(
    metabolite_group: ty.MetaboliteGroup | list[str],
    hue: ty.FilterColor,
    start_diagnosis: ty.Diagnosis,
    sample_source: ty.Sample = ty.Sample.Gut,
    end_diagnosis: ty.Diagnosis = ty.Diagnosis.All,
    biological_sex: ty.BioSex = ty.BioSex.All,
) -> pd.DataFrame:
    assert start_diagnosis != ty.Diagnosis.All

    data, hue_value, variables, _ = _prepare_data_and_variables(
        metabolite_group=metabolite_group,
        hue=hue,
        start_diagnosis=start_diagnosis,
        sample_source=sample_source,
        end_diagnosis=end_diagnosis,
        biological_sex=biological_sex,
    )

    results = []

    for variable, title in list(variables.items())[2:]:
        model_formula = f"{variable} ~ C({hue_value})"
        model = ols(model_formula, data=data).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        p_value = round(anova_table.loc[f"C({hue_value})", "PR(>F)"], 4)
        ss_effect = anova_table.loc[f"C({hue_value})", "sum_sq"]
        ss_total = anova_table["sum_sq"].sum()
        eta_squared = ss_effect / ss_total  # Effect size η²

        if p_value < 0.05:
            group_names = data[hue_value].unique()
            pairs = list(combinations(group_names, 2))
            pairwise_pvalues = []
            pairwise_effect_sizes = []
            group_means = {
                group: data[data[hue_value] == group][variable].mean()
                for group in group_names
            }

            for group1, group2 in pairs:
                group1_data = data[data[hue_value] == group1][variable]
                group2_data = data[data[hue_value] == group2][variable]
                _, p = f_oneway(group1_data, group2_data)
                pairwise_pvalues.append(p)

                # Calculate Cohen's d
                mean_diff = np.mean(group1_data) - np.mean(group2_data)
                pooled_std = np.sqrt(
                    (
                        (len(group1_data) - 1) * np.var(group1_data, ddof=1)
                        + (len(group2_data) - 1) * np.var(group2_data, ddof=1)
                    )
                    / (len(group1_data) + len(group2_data) - 2)
                )
                cohens_d = mean_diff / pooled_std
                pairwise_effect_sizes.append(cohens_d)

            # Bonferroni correction
            _, corrected_pvalues, _, _ = multipletests(
                pairwise_pvalues, method="bonferroni"
            )
            significance_level = 0.05 / len(pairwise_pvalues)

            # Store results
            for (group1, group2), corrected_p, cohens_d in zip(
                pairs, corrected_pvalues, pairwise_effect_sizes
            ):
                if corrected_p < significance_level:
                    higher_mean_group = (
                        group1 if group_means[group1] > group_means[group2] else group2
                    )
                    results.append(
                        {
                            "Variable": title,
                            "Comparison": f"{group1} vs {group2}",
                            "Corrected p-value": round(corrected_p, 4),
                            "η²": round(eta_squared, 4),
                            "Cohen's d": round(cohens_d, 4),
                            "Higher Mean Group": higher_mean_group,
                        }
                    )
        if p_value > 0.05:
            results.append(
                {
                    "Variable": title,
                    "Corrected p-value": p_value,
                }
            )

    results_df = pd.DataFrame(results)
    return results_df


def plot_age_vs_variable(
    metabolite_group: ty.MetaboliteGroup | list[str],
    hue: ty.FilterColor,
    start_diagnosis: ty.Diagnosis,
    sample_source: ty.Sample,
    end_diagnosis: ty.Diagnosis = ty.Diagnosis.All,
    biological_sex: ty.BioSex = ty.BioSex.All,
    color_config: ty.FilterColor = ty.FilterColor.Sex,
):
    """
    Plots scatter plots of "AGE" against each variable with regression lines.

    Arguments:
    - metabolite_group: This takes on a class or a list of column names from the original dataset.
    - hue: The column to group data by.
    - start_diagnosis: The starting diagnosis filter.
    - end_diagnosis: The ending diagnosis filter (default is Diagnosis.All).
    - biological_sex: The gender filter (default is BioSex.All).
    - color_config: FilterColor to determine plot colors.
    """
    assert start_diagnosis != ty.Diagnosis.All

    data, hue_value, variables, hist_color = _prepare_data_and_variables(
        metabolite_group=metabolite_group,
        hue=hue,
        start_diagnosis=start_diagnosis,
        sample_source=sample_source,
        end_diagnosis=end_diagnosis,
        biological_sex=biological_sex,
        color_config=color_config,
    )

    columns = 4
    rows = math.ceil(len(variables) / columns)  # Determine the required number of rows
    fig, axes = plt.subplots(nrows=rows, ncols=columns, figsize=(6 * columns, 6 * rows))
    axes = axes.ravel()

    for i, (variable, title) in enumerate(list(variables.items())[2:]):
        sns.scatterplot(
            data=data,
            x="AGE",
            y=variable,
            hue=hue_value,
            ax=axes[i],
            palette=hist_color,
        )
        axes[i].set_title(title)
        axes[i].set_xlabel("Age (years)")
        axes[i].set_ylabel("Metabolite Concentration")

    # Hide any unused subplots
    for j in range(i + 1, rows * columns):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.show()


def _prepare_data_and_variables(
    metabolite_group: ty.MetaboliteGroup | list[str],
    hue: ty.FilterColor,
    start_diagnosis: ty.Diagnosis,
    sample_source: ty.Sample,
    end_diagnosis: ty.Diagnosis = ty.Diagnosis.All,
    biological_sex: ty.BioSex = ty.BioSex.All,
    color_config: ty.FilterColor = ty.FilterColor.Sex,
):
    """
    Prepares the parameters for the main functions.
    """
    data = fu.add_enddx_total_to_baseline(
        metabolites="gut_metabolites_result.csv",
        diagnoses="All_Subjects_DXSUM_08Oct2024.csv",
        gender="ADNIMERGE_08Oct2024.csv",
        sample_source=sample_source,
        metabolite_group=metabolite_group,
        start_diagnosis=start_diagnosis,
        end_diagnosis=end_diagnosis,
    )

    match biological_sex:
        case ty.BioSex.Female:
            data = data.loc[data["PTGENDER"] == "Female"]
        case ty.BioSex.Male:
            data = data.loc[data["PTGENDER"] == "Male"]

    match hue:
        case ty.FilterColor.Sex:
            hue_value = "PTGENDER"
        case ty.FilterColor.Diagnosis:
            hue_value = "END_DIAGNOSIS"

    variables = fu._name_variables(
        metabolite_group=metabolite_group, sample_source=sample_source
    )

    hist_color = fu.get_color_config(color_config=color_config)

    return data, hue_value, variables, hist_color
