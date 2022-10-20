Data Acquisition and Preprocessing
==================================

Set up Computing Environment
----------------------------

**Note, this code was tested on Ubuntu Linux.** These instructions are
intended to minimize setup time and (if necessary) costs associated with
hourly compute resources (~$1/hour on
[DigitalOcean](https://www.digitalocean.com/pricing/calculator/) for 128
GB RAM, 12 core processor.)

**Install GDAL**

> `sudo apt update` `sudo apt upgrade`

Install dependencies:
&gt;`sudo apt-get install gdal-bin`

**Clone the repository (from the root directory)**

> `git clone https://github.com/dankovacek/hysets_validation`

Change directories to the `hysets_validation` folder:
&gt;`cd hysets_basin_delineation`

### Install Python package manager (pip)

If not automatically installed, install Python and virtualenv (assuming
Python3 is installed by default on a linux distribution):

> `sudo apt install python3.8-venv pip`

Create Python 3.8+ virtual environment at the root level directory:

> `python3 -m venv env/`

Activate the virual environment:  
&gt;`source env/bin/activate`

Install Python packages:  
&gt;`pip install -r requirements.txt`

Data Acquisition and Processing
-------------------------------

### EarthEnv DEM90 (Robinson, Regetz, and Guralnick 2014)

HYSETS used EarthEnv DEM90 (~90m resolutoin) data to derive basins and
physical / topographical attributes. The tiles are available at
[earthenv.org/DEM](https://earthenv.org/DEM). We can download the tiles
and merge them into a virtual raster with gdal using the following
script in `setup_scripts`:

> `python get_EarthEnv_DEM.py`

Links to invidivual DEM tiles look like the following:
&gt;`http://mirrors.iplantcollaborative.org/earthenv_dem_data/EarthEnv-DEM90/EarthEnv-DEM90_N55W125.tar.gz`

The resulting .vrt mosaic should look like below (red outline added to
emphasize land mass bounds.):

![DEM Mosaic of BC and administrative boundary
regions](../images/CANUSA_Mex_bounds.png)

### HYSETS

The complete set of files associated with the HYSETS paper is ~14.6 GB
and can be accessed [here](https://osf.io/rpc3w/). However, for
validation we only need the results file, basin geometry, and station
locations. The results file is included in the repo, and we can download
the rest with the following script, found in `setup_scripts/`:

> `cd setup_scripts` `python get_HYSETS_data.py`

The basin polygon file (~15GB) could take 15-20 minutes on it’s own to
download, so open a second terminal window and continue the setup
process while this file is downloading.

The process of delineating basins requires the specification of a pour
point for each station. The accuracy of basin delineation is highly
sensitive to the specification of pour points. Pour points are not
specified for the polygons derived for this study, so these must be
derived from WSC and USGS data. Pour points for WSC stations were
recently published and these are downloaded in the next step. You may
find that station locations do not yield accurate catchment delineations
when using station locations.

Next we will download basin polygons from the WSC. These files contain
both station locations and pour points which produce much better results
when validating the catchment polygon.

<!-- ### WSC Hydrometric Station Catchment Polygons and Metadata

An updated set of basin polygons from WSC, published in December 2021, can be downloaded and processed by the following script located in `setup_scripts`.  You can customize the region of interest by modifying the list of basin prefixes, i.e. BC is covered by 07, 08, 09, and 10.

>`python get_WSC_data.py`

The `source_data` folder contains the WSC station metadata file of active and historic hydrometric stations `WSC_Stations_2020.csv`.  


Next we choose one or both DEM sources to download tiles and construct a DEM mosaic file covering the study area. -->

Basin Delineation and Polygon Validation
----------------------------------------

This step represents the heavy lifting where large regions of DEM such
as the Liard, Peace, and and Fraser River basins are processed into flow
direction and flow accumulation at the highest available resolution.
Three tools were evaluated for the step of validating HYSETS basins:
[Whiteboxtools](python%20process_dem_by_basin.py),
[RichDEM](https://richdem.readthedocs.io/en/latest/), and
[Pysheds](https://mattbartos.com/pysheds/). Each have distinct feature
sets, but Pysheds was used here for the step of delineating a large set
of basins.

<!-- Note: the breach [depression function](https://jblindsay.github.io/ghrg/Whitebox/Help/BreachDepressions.html) run on the DEM is a bottleneck step.   -->

The final step is to validate the basin attributes derived in HYSETS (or
other dataset) using the set of stations whose catchment boundaries
intersect BC. The manual basin delineation step is the most
computationally intensive step of the validation process, and it’s
executed with the script `setup_scripts/delineate_basins.py`

> `setup_scripts/python pysheds_derive_basin_polygons.py`

Additional Notes
----------------

To automate citation formatting for the README document.

> `pandoc -t markdown_strict -citeproc README-draft.md -o README.md --bibliography bib/bibliography.bib`

Robinson, Natalie, James Regetz, and Robert P Guralnick. 2014.
“EarthEnv-Dem90: A Nearly-Global, Void-Free, Multi-Scale Smoothed, 90m
Digital Elevation Model from Fused Aster and Srtm Data.” *ISPRS Journal
of Photogrammetry and Remote Sensing* 87: 57–67.
