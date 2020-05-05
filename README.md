# Retrieve, format, and impute Tara Oceans environmental data 
Combines climatological means of iron, phosphate, nitrate from the [MIT DARWIN model](https://darwinproject.mit.edu/) accessed from [Simons Collaborative Marine Atlas Project](https://simonscmap.com/), supplementary data from Sunagawa supplementary material found [here](http://ocean-microbiome.embl.de/companion.html), and environmental measurements (hydrography, biogeochemistry) from [PANGAEA](https://doi.pangaea.de/10.1594/PANGAEA.875576)

The raw data is not provided here so as to not distribute anything outside of license permissions. To run this pipeline you will need the Sunagawa [supplementary excel table](http://ocean-microbiome.embl.de/data/OM.CompanionTables.xlsx) as well as the complete data sensors and summarized [PANGAEA dataset from here](https://doi.pangaea.de/10.1594/PANGAEA.875576) (scroll to bottom of page, click "Download dataset as tab-delimited text").

You will also need to [register with Simons CMAP](https://simonscmap.com/register), obtain and API key, then download the `cmap4r` [package](https://simonscmap.github.io/cmap4r/index.html) in order to access data from the MIT DARWIN model.

You can download the Longhurst Codes and identifiers for Tara samples (or run it yourself) from [here](https://github.com/slhogle/Longhurst-Province-Finder)

Relative ecotype abundances of SAR11 and Prochlorococcus are available here in the data repository. This data was originally published in 
Becker et al. 2019.

## ATTRIBUTION
Please only download and use data if you agree to respective license terms for each data source!

If you use any of the data mentioned or described in this repository it is ESSENTIAL!! that you properly cite the data sources. These include (but are not limited to):

### Longhurst codes, identifiers, and descriptions
Longhurst, A. R. 1998. Ecological Geography of the Sea, Academic Press.

https://github.com/slhogle/Longhurst-Province-Finder

### Prochlorococcus ecotype abundances from metagenomes
Becker, J. W., S. L. Hogle, K. Rosendo, and S. W. Chisholm. 2019. Co-culture and biogeography of Prochlorococcus and SAR11. ISME J. doi:10.1038/s41396-019-0365-4

### MIT DARWIN Model output
Simons Collaborative Marine Atlas project (CMAP) https://simonscmap.com/

Dutkiewicz, S., A.E. Hickman, O. Jahn, W.W. Gregg, C.B. Mouw, and M.J. Follows, 2015:  Capturing optically important constituents and properties in a marine biogeochemical and ecosystem model. Biogeoscience, 12, 4447-4481 doi:10.5194/bg-12-4447-2015, https://doi.org/10.5194/bg-12-4447-2015

Forget, G., Campin, J.-M., Heimbach, P., Hill, C. N., Ponte, R. M., and Wunsch, C.: ECCO version 4: an integrated framework for non-linear inverse modeling and global ocean state estimation, Geosci. Model Dev., 8, 3071-3104, https://doi.org/10.5194/gmd-8-3071-2015, 2015

Forget, G., D. Ferreira, and X. Liang, 2015: On the observability of turbulent transport rates by argo: supporting evidence from an inversion experiment. Ocean Science, 11, 839–853, doi:10.5194/os-11-839-2015

Forget, G. and R. Ponte, 2015: The partition of regional sea level variability. Progress in Oceanography, 137, 173–195, https://doi.org/10.1016/j.pocean.2015.06.002

Forget, G., 2018: Initial, preliminary version of the CBIOMES-global model setup and documentation (Version v0.0.1). Zenodo. http://doi.org/10.5281/zenodo.1343303

Forget, G., 2019: Update MITgcm & DarwinProject elements (Version v0.1.0). Zenodo. http://doi.org/10.5281/zenodo.2653669
Ward, B.A., S. Dutkiewicz, O. Jahn, and M.J. Follows, 2012: A size-structured food-web model for the global ocean. Limnol. Oceanogr., 57, 1877-1891. https://aslopubs.onlinelibrary.wiley.com/doi/abs/10.4319/lo.2012.57.6.1877

### Tara Oceans datasets
Sunagawa, S., L. P. Coelho, S. Chaffron, and others. 2015. Structure and function of the global ocean microbiome. Science 348: 1261359.

Guidi, Lionel; Picheral, Marc; Pesant, Stephane; Tara Oceans Consortium, Coordinators; Tara Oceans Expedition, Participants (2017): Environmental context of all samples from the Tara Oceans Expedition (2009-2013), about sensor data in the targeted environmental feature. PANGAEA, https://doi.org/10.1594/PANGAEA.875576

### Code
If any of this code is helpful for you please also consider citing this repository. [![DOI](https://zenodo.org/badge/261384041.svg)](https://zenodo.org/badge/latestdoi/261384041)

## WORKFLOW
[Rendered markdown available here](bin/formatting_pipeline.md)