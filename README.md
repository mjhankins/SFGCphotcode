# SFGCphotcode
* Python code developed for constructing the source catalog for the SOFIA/FORCAST Galactic Center Legacy Program. 

The code here includes steps for source detection and photometry on a field by field basis. There are also notebooks for combining the results into the final catalog.

Typical work flow for using this code is as follows:

1. Specify observations in 'config.py' file. An example configuration file is provided named 'config_example.py' but will need to be renamed and contents modified to run on a new machine. 
2. Run SourceCatalog_Detect.ipynb
3. Run SourceCatalog_ApPhot.ipynb
4. Run SourceCatalog_CombineFields.ipynb
5. Run SourceCatalog_SSTGC_Xmatch.ipynb

Note that steps 2 and 3 have python script versions with the same name. These scripts will perform batch runs for all files specified in the configuration file. 

The final catalog data is available in the 'Data' directory in this repository. This includes a version of the catalog with image cutouts for all sources. 

The catalog paper is currently submitted to ApJ, and upon publication, we will update this page to include links to the paper which provides further description of the files here. 

This code can be referenced using the following information:

[![DOI](https://zenodo.org/badge/376088869.svg)](https://zenodo.org/doi/10.5281/zenodo.11459087)



