# **Hotspot propensity across mutational processes**
This repository contains the source code to reproduce hotspot propensity figures from the manuscript: 

*Arnedo-Pac C, Muiños F, Gonzalez-Perez A, Lopez-Bigas N. Hotspot propensity across mutational processes. bioRxiv 2022.09.14.507952; doi: [10.1101/2022.09.14.507952](https://www.biorxiv.org/content/10.1101/2022.09.14.507952v2)*

## Content
The repository contains the data and code to reproduce the main and Extended View figures figures of the manuscript. You can also check the resulting PNG images from running the code. 
- [Figure 1](https://github.com/bbglab/hotspot_propensity/tree/main/figures/main_figures/figure_1)
- [Figure 2](https://github.com/bbglab/hotspot_propensity/tree/main/figures/main_figures/figure_2)
- [Figure 3](https://github.com/bbglab/hotspot_propensity/tree/main/figures/main_figures/figure_3)
- [Figure 4](https://github.com/bbglab/hotspot_propensity/tree/main/figures/main_figures/figure_4)
- [Figure 5](https://github.com/bbglab/hotspot_propensity/tree/main/figures/main_figures/figure_5)
- [Extended View Figure 1](https://github.com/bbglab/hotspot_propensity/tree/main/figures/extended_view_figures/figure_EV1)
- [Extended View Figure 2](https://github.com/bbglab/hotspot_propensity/tree/main/figures/extended_view_figures/figure_EV2)

## How to run
#### Requirements

This code has been developed in Python 3.6 and Jupyter Notebook version 5.0.0. You can check the list of additional requirements at [requirements.txt](https://github.com/bbglab/hotspot_propensity/blob/main/requirements.txt) file. 

#### Step 1. Create a conda environment and install the requirements
You can directly create a conda environment as: 
```sh
~$ conda create --name hotspot_propensity_env --file requirements.txt
```

#### Step 2. Get a copy of hotspot_propensity repository
For example, you can clone the repository as:
```sh
~$ git clone git@github.com:bbglab/hotspot_propensity.git
```
#### Step 3. Download and uncompress source data
The source data is available within this [Github](https://github.com/bbglab/hotspot_propensity) and a Zenodo [10.5281/zenodo.10004773](https://doi.org/10.5281/zenodo.10004773) repositories.

1. Download hotspots identified across cancer types from Zenodo [10.5281/zenodo.10004773](https://doi.org/10.5281/zenodo.10004773). 
2. Move `hotspots.zip` file to `./data` directory and uncompress. 
3. The resulting `./data` directory should look like this: 

```sh
~$ ls ./data -1
colonic_crypts_mutations
ctcf
EV_datasets
expected_hotspot_propensity
fraction_samples_with_hotspot.json
genomic_bin_data
germline_mutations
hotspot_propensity_1000iter_100samples_100-300muts.txt.gz
hotspots
hotspots_SBS_prob__COADREAD.txt.gz
hotspots_SBS_prob__ESOPHA_STOMACH.txt.gz
hotspots_SBS_prob__SKCM.txt.gz
methylation
sample_signatures_fraction_activity_inside_outside.json
signatures_fold_change_inside_outside.txt
total_hotspots_per_sample.json
total_mutations_per_sample.json
```
#### Step 4. Activate the conda environment and run Jupyter Notebook
```sh
~$ conda activate hotspot_propensity_env
~$ jupyter notebook
```
Run the notebook of your interest. Figures will appear in the running directory of the notebook. 

## Additional information
You can read all details of the manuscript at: Arnedo-Pac C, Muiños F, Gonzalez-Perez A, Lopez-Bigas N. Hotspot propensity across mutational processes. bioRxiv 2022.09.14.507952; doi: [10.1101/2022.09.14.507952](https://www.biorxiv.org/content/10.1101/2022.09.14.507952v2)

HotspotFinder algorithm is described in Materials and Methods section and Appendix Note 1 of the manuscript and its repository. You can download HotspotFinder at [bitbucket.org/bbglab/hotspotfinder](https://bitbucket.org/bbglab/hotspotfinder/src/master/)

## How to cite
Arnedo-Pac C, Muiños F, Gonzalez-Perez A, Lopez-Bigas N. Hotspot propensity across mutational processes. bioRxiv 2022.09.14.507952; doi: [10.1101/2022.09.14.507952](https://www.biorxiv.org/content/10.1101/2022.09.14.507952v2)

## License
This code is available to the general public subject to certain conditions described in its [LICENSE](https://github.com/bbglab/hotspot_propensity/blob/main/LICENSE). 


