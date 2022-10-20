BC Streamflow Monitoring Network Analysis
=========================================

Introduction
------------

The importance of validating inputs to any modelling effort is
important. Large sample hydrometric datasets have a key role as input in
hydrological studies (Addor et al. 2020; Gupta et al. 2014; Klingler,
Schulz, and Herrnegger 2021) in both process-based and machine learning
approaches. Kratzert et al. (2019) used the large sample hydrologic
dataset CAMELS (Addor et al. 2017) in a machine learning model for
rainfall-runoff prediction in ungauged basins. Other efforts in
generating large sample datasets attribute sets similar to CAMELS
quickly followed (Alvarez-Garreton et al. 2018; Coxon et al. 2020;
Fowler et al. 2021).

Arsenault et al. (2020) developed [HYSETS](https://osf.io/rpc3w/), a
large-sample dataset of over 14K hydrometric monitoring stations in
North America and Mexico. This repo presents a method of replicating the
attributes derived for the Arsenault et al. (2020) study, from
collecting the source spatial data, to deriving catchments, to deriving
basin attributes and comparing to the values derived in the study.

Basin delineation is critical to acquiring representative values for
many attribute indices, and here we discuss several challenges in
automated basin delineation and highlight the information added to the
process manually by experienced GIS practitioners, geographers, and
hydrologists.

Notes
-----

Setup
-----

See the [README.md under
`setup_scripts/`](https://github.com/dankovacek/hysets_basin_delineation/tree/main/setup_scripts)
for setup of the validation scripting.

References
----------

Addor, Nans, Hong X Do, Camila Alvarez-Garreton, Gemma Coxon, Keirnan
Fowler, and Pablo A Mendoza. 2020. “Large-Sample Hydrology: Recent
Progress, Guidelines for New Datasets and Grand Challenges.”
*Hydrological Sciences Journal* 65 (5): 712–25.

Addor, Nans, Andrew J Newman, Naoki Mizukami, and Martyn P Clark. 2017.
“The Camels Data Set: Catchment Attributes and Meteorology for
Large-Sample Studies.” *Hydrology and Earth System Sciences* 21 (10):
5293–5313.

Alvarez-Garreton, Camila, Pablo A Mendoza, Juan Pablo Boisier, Nans
Addor, Mauricio Galleguillos, Mauricio Zambrano-Bigiarini, Antonio Lara,
et al. 2018. “The Camels-Cl Dataset: Catchment Attributes and
Meteorology for Large Sample Studies–Chile Dataset.” *Hydrology and
Earth System Sciences* 22 (11): 5817–46.

Arsenault, Richard, François Brissette, Jean-Luc Martel, Magali Troin,
Guillaume Lévesque, Jonathan Davidson-Chaput, Mariana Castañeda
Gonzalez, Ali Ameli, and Annie Poulin. 2020. “A Comprehensive,
Multisource Database for Hydrometeorological Modeling of 14,425 North
American Watersheds.” *Scientific Data* 7 (1): 1–12.

Barnes, Richard. 2016. “Parallel Priority-Flood Depression Filling for
Trillion Cell Digital Elevation Models on Desktops or Clusters.”
*Computers & Geosciences* 96: 56–68.

Coxon, Gemma, Nans Addor, John P Bloomfield, Jim Freer, Matt Fry, Jamie
Hannaford, Nicholas JK Howden, et al. 2020. “CAMELS-Gb:
Hydrometeorological Time Series and Landscape Attributes for 671
Catchments in Great Britain.” *Earth System Science Data* 12 (4):
2459–83.

Fowler, Keirnan JA, Suwash Chandra Acharya, Nans Addor, Chihchung Chou,
and Murray C Peel. 2021. “CAMELS-Aus: Hydrometeorological Time Series
and Landscape Attributes for 222 Catchments in Australia.” *Earth System
Science Data* 13 (8): 3847–67.

Gleeson, Tom, Nils Moosdorf, Jens Hartmann, and LPH Van Beek. 2014. “A
Glimpse Beneath Earth’s Surface: GLobal Hydrogeology Maps (Glhymps) of
Permeability and Porosity.” *Geophysical Research Letters* 41 (11):
3891–8.

Gupta, Hoshin Vijai, C Perrin, G Blöschl, A Montanari, R Kumar, M Clark,
and Vazken Andréassian. 2014. “Large-Sample Hydrology: A Need to Balance
Depth with Breadth.” *Hydrology and Earth System Sciences* 18 (2):
463–77.

Huscroft, Jordan, Tom Gleeson, Jens Hartmann, and Janine Börker. 2018.
“Compiling and Mapping Global Permeability of the Unconsolidated and
Consolidated Earth: GLobal Hydrogeology Maps 2.0 (Glhymps 2.0).”
*Geophysical Research Letters* 45 (4): 1897–1904.

Klingler, Christoph, Karsten Schulz, and Mathew Herrnegger. 2021.
“LamaH-Ce: LArge-Sample Data for Hydrology and Environmental Sciences
for Central Europe.” *Earth System Science Data* 13 (9): 4529–65.

Kratzert, Frederik, Daniel Klotz, Mathew Herrnegger, Alden K Sampson,
Sepp Hochreiter, and Grey S Nearing. 2019. “Toward Improved Predictions
in Ungauged Basins: Exploiting the Power of Machine Learning.” *Water
Resources Research* 55 (12): 11344–54.

Robinson, Natalie, James Regetz, and Robert P Guralnick. 2014.
“EarthEnv-Dem90: A Nearly-Global, Void-Free, Multi-Scale Smoothed, 90m
Digital Elevation Model from Fused Aster and Srtm Data.” *ISPRS Journal
of Photogrammetry and Remote Sensing* 87: 57–67.
