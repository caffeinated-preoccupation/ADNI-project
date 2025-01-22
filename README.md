# Gut and Serum Metabolites and Alzheimer’s: Exploring the Link Between Metabolite Levels and the Risk of Developing Alzheimer’s

This project aims to explore the link between the future Alzheimer's risk of a patient given some metabolite metrics. The data used for this compilation of exploratory data analysis was downloaded from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu).

## Installation

To make sure all the package dependencies work, use pip install to install the requirements as listed in the requirements file.

```bash
pip install -r requirements.txt
```

## Usage

There are different python files contained in this project and the analysis of corresponding metabolite levels require the import of each file. Below are the descriptions of each of the relevant python files:

### gut_metabolites_labels.py

```python
import gut_metabolites_labels
```

Contains all the dictionaries used for gut metabolites; includes the mapping of the column names in the data to the actual column names, as indicated in the data dictionary provided by ADNI. This file was imported into the shared_functions.py.

### serum_metabolites_labels.py

```python
import serum_metabolites_labels
```

Contains all the dictionaries used for serum metabolites; includes the mapping of the column names in the data to the actual column names, as indicated in the data dictionary provided by ADNI. This file was imported into the shared_functions.py.

### shared_types.py
```python
import shared_types.py
```

Contains all Enum classes and frozen dataclass used in the functions within this project.

### shared_functions.py
```python
import shared_functions.py
```
Contains all the functions that were used for the analysis of both serum and gut metabolites. 

#### Available Functions
- subset_start_end_diagnosis: returns a subset of the data based on start and end diagnosis filter 
- add_enddx_total_to_baseline: adds the end diagnosis and the sum of the metabolites as columns to the data 
- display_pairwise_results_on_plot: display the pairwise statistic on subplots 
- get_color_config: returns the color palette to be used for the subplots 
- plot_sunburst_dist: returns a sunburst plot based on selected data path 

### gut_metabolites.py

```python
import gut_metabolites as gm
```
Contains all the dictionaries used for the analysis of gut metabolites.

#### Available Functions
- prepare_data_and_variables: returns a tuple of parameters 
- plot_hist_kde: plots the histogram and density plot with corresponding statistic results 
- display_effect_sizes_and_cohens: returns a df of the statistic results for effect size 
- plot_age_vs_variable: returns multiple subplots of age with corresponding variables 

### serum_metabolites.py

```python
import serum_metabolites as sm
```
Contains all the dictionaries used for the analysis of serum metabolites.

#### Available Functions
- prepare_data_and_variables: returns a tuple of parameters 
- plot_hist_kde: plots the histogram and density plot with corresponding statistic results 
- display_effect_sizes_and_cohens: returns a df of the statistic results for effect size 
- plot_age_vs_variable: returns multiple subplots of age with corresponding variables 


## Examples

The accompanying Jupyter notebooks show examples on how to use the data analysis pipelines. These notebooks do not involve all the codes used for this project as they are repeated.

### gut_metabolites.ipynb
Contains examples for the analysis of the gut metabolites.

### serum_metabolites.ipynb
Contains examples for the analysis of the serum metabolites.

## Contributing
As this is a study project that is completed, the project is not open for contribution.  

## License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)