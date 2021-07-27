<img src="https://github.com/iquasere/MOSCA/blob/master/mosca_logo.png" align="center" height="300">

Logo by [Sérgio A. Silva](https://www.ceb.uminho.pt/People/Details/64888072-5cde-42da-b7e5-691d380cefb2)

# Meta-Omics Software for Community Analysis (MOSCA)

**MOSCA** (portuguese for fly) is a pipeline designed for performing metagenomics (MG) and metatranscriptomics (MT) integrated data analyses, in a mostly local and fully automated workflow.

[![CI](https://github.com/iquasere/MOSCA/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/iquasere/MOSCA/actions/workflows/main.yml)

## MOSCA's wiki

Information about MOSCA is primarily available through the [wiki](https://github.com/iquasere/MOSCA/wiki), where detailed information about the pipeline can be found.

It includes the following sections:
* [Installing-and-running-MOSCA](https://github.com/iquasere/MOSCA/wiki/Installing-and-running-MOSCA)
* [Parameters of MOSCA](https://github.com/iquasere/MOSCA/wiki/Parameters-of-MOSCA)
* [Partial runs](https://github.com/iquasere/MOSCA/wiki/Partial-runs)
* [Project and community](https://github.com/iquasere/MOSCA/wiki/Project-and-community)
* [Technical documentation](https://github.com/iquasere/MOSCA/wiki/Technical-documentation)

## Quick setup of MOSCA

MOSCA can be quickly setup by downloading it and compiling from source code with
```
git clone https://github.com/iquasere/MOSCA.git
bash MOSCA/workflow/envs/install.bash
conda activate mosca
```
Two configuration files are required for MOSCA, which can be obtained with MOSCA's GUI: [MOSGUITO](https://iquasere.github.io/MOSGUITO/). After obtaining them, MOSCA can be run with this command
```
mosca.py -c config.json
```
where ```config.json``` is the configuration file obtained with MOSGUITO. Do note that the ```experiments``` parameter should indicate the path to the ```experiments``` file, the second configuration file for MOSCA.

If the wiki is not enough to elucidate about MOSCA, or if you find any problem related to using MOSCA, don't hesitate in asking, either through an [issue](https://github.com/iquasere/MOSCA/issues) or by directly contacting the developer of MOSCA at jsequeira@ceb.uminho.pt.

If you would like to contribute to this project, a [pull request](https://github.com/iquasere/MOSCA/pulls) is the way to go!
