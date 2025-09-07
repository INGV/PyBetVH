# PyBetVH

PyBetVH is an open and cross-platform implementation of BET_VH model. The main purpose of this software is to provide a graphically supported computation of long-term probabilities of volcanic hazardous phenomena (i.e., lava flows, tephra fall, pyroclastic flows, lahars, etc.) through the Bayesian Event Tree model for Volcanic Hazard (BET_VH, Marzocchi et al., 2010). The model represents a flexible tool to provide probabilities of any specific event at which we are interested in, by merging all the available information, such as models, a priori beliefs, and past data. It is mainly based on a Bayesian procedure, in order to quantify the aleatory and epistemic uncertainty characterizing the impact of volcanic eruptions in terms of hazard assessment. The method deals only with long-term forecasting, therefore in principle it can be useful for land use planning. This version introduces as main outputs the Hazard Curves (HC), showing the exceeding probability as function of the chosen intensity measure (e.g., the tephra load in kg/m2). The HC are provided together with their uncertainties (percentiles) through Bayesian inference and the user can explore the hazard results visualizing both Probability and Hazard Maps. 

More datails on the PyBetVH software can be found in Tonini et al. (2015). A general user guide for all PyBet tools can be found at [https://theghub.org/wiki/PyBetToolsUserGuide](https://theghub.org/wiki/PyBetToolsUserGuide) (formerly Vhub cyberinfrastructure platform).

**References:** 
 - *Marzocchi W., Sandri L., Selva J. (2010) BETVH: a probabilistic tool for long-term volcanic hazard assessment, Bull. Volcanol., 72, 705-716, [doi:10.1007/s00445-010-0357-8](https://link.springer.com/article/10.1007%2Fs00445-010-0357-8)* 
 - *Tonini R., Sandri L., Thompson M. A. (2015) PyBetVH: a Python tool for probabilistic volcanic hazard assessment and for generation of Bayesian hazard curves and maps, Comput. Geosci., [doi:10.1016/j.cageo.2015.02.017](https://www.sciencedirect.com/science/article/pii/S0098300415000515)*


## Requirements
This version of PyBetVH requires Python 3.x and the following Python modules/libraries:
 - wxpython
 - numpy
 - scipy
 - matplotlib
 - pillow (formerly known as PIL)


## Installation
Download this repository by clicking the button placed on the top-right of this page, under the project title (the one with the small cloud icon with a downward arrow) and unzip the archive wherever you prefer in your local computer (referred as `/path_to_pybetvh/` from hereafter).

Alternatively, you can make a clone of the project:
```
git clone git@github.com:INGV/PyBetVH.git
```

Now, in order to run the tool, you need a Python interpreter and the packages listed above (compatible with your Python interpreter version).
Different approaches can be adopted to satisfy these requierments, depending on the experience of the user with Python programming and on the used operating system.

NOTE: The software is developed and tested mainly on Debian-based Linux distribution, meaning that software issues due to different operating systems could be missed by the developers.


#### Using a Conda environment
A suggested (and cross-platform) way to run PyBetUnrest is to use a [Conda](https://conda.io/en/latest/) environment, in particular [Miniconda](https://conda.io/en/latest/miniconda.html). This will allow to build a specific Python environment separated from your system libraries.
Moreover, Conda is a package manager, so it will take care for you of check the correct dependencies among packages.

NOTE: Users that are working with the much more complete [Anaconda distribution](https://www.anaconda.com/), should be able to run the PyBetVH tool without installing Miniconda. 

After having downloaded the last Miniconda (if you do not have particoular needs, the Miniconda Python3 version is suggested) for your own platform from [here](https://conda.io/en/latest/miniconda.html), install it, by accepting the licence and following the default options suggested by the installer. The installer will create a miniconda3/ (depending on the installed version of Miniconda) folder in a given default path (hereafter, I'll refer to it with `/path_to_miniconda/`), depending on your operating system.
Linux (and OSX) users can run the following command lines to download and install last Miniconda version: 

```
wget -c http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

 - **Linux and Mac/OSX users**
The installer should add the Miniconda environment to your PATH variable. 
You can check it opening the terminal and running `echo $PATH`
If Miniconda is not in your PATH variable, you can temporary add it by running:
`export PATH="path_to_miniconda/bin:$PATH" `
or adding the previous line to the .profile or .bash_profile.

 - **Windows users**
In Windows the installer suggests to not add the PATH since they can use directly the Anaconda Prompt from the WIndows munu, which load the Conda environment. 


Then you have to create the specific environment to run PyBetVH. To do this, open a terminal (Windows users have to do it from the Miniconda Prompt) and type:
```
conda env create -f /path_to_pybetvh/pybet.yml
```

If this command returns an error, you can manually create the environment as follows (here python 3.9 is used, but any 3.x version considered `stable` should work):
```
conda create --name pybet python=3.9
conda activate pybet
conda install numpy pillow matplotlib wxpython
```

If `conda activate pybet` returns an error, try `source activate pybet`, depending on versions of conda itself.

NOTE: if you find some errors during the installation of wxpython library, try adding the conda-forge channel for conda by typing this command out of the _pybet_ environment:

```
conda config --add channels conda-forge
```


Once the installation worked fine, the user can activate the pybet environment. 
When the environment is activate, the prompt displays its name between brackets:
```
(pybet) prompt: 
```

The shell will continue to work regularly but using the Python interpreter and the libraries of the activated environment. 
To deactivate it just type:
 - `conda deactivate`
or
 - `source deactivate`
or close the shell and open another one.
  
Once you have the pybet environment activated, you can run the tool:
```
(pybet) prompt: python /path_to_pybetvh/src/betvh.py
```
