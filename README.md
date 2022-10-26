# Analysis tools for Minflux microscopy data

Tools are available in Matlab and Python (both languages
may not be equally supported). Scripts are under the MIT
license (see LICENSE).

## Import of Abberior Instruments Minflux microscopy data

This section goes into import of AI Minflux microscopy
(https://abberior-instruments.com/products/minflux/) data.

### Exported localizations format (May 2022)

The exported Minflux localization data format was subject
to some changes between 2021 and 2022. A Matlab (*.mat),
a Python (*.npy) or a JSON (*.json) export option exists.

They all represent localization as nested arrays of structs.

The top level consists of fields:

- gri  ?
- itr  (itself a struct see below)
- sqi  ?
- tid  trace id (unique for consecutive localizations of
  the same mol.)
- tim  time of localization (s) of last iteration (equivalent
  to itr.tic(:, end) / 4e-7)
- vld  true if regarded as valid, otherwise false

All fields except for itr are 1xN vectors for N primary
localization events. (Software allows to include, exclude
invalid localization attempts.)

Field itr is a struct with fields:

- cfr  center frequency ratio (equals efc / efo)
- dcr  detector channel ratio (ratio of Cy5 near to Cyr5 far
  detector counts - or the other way around)
- dmz  deformable mirror z-position
- eco  effective counts periphery (circumference of beampattern)
- ecc  effective counts center (middle position of beampattern)
- efo  like eco but frequency (Hz)
- efc  like ecc but frequency (Hz)
- ext  beampattern diameter (scales with wavelength and with
  L in sequence)
- eox  EOD position relative to Galvo center in x
- eoy  .. in y
- fbg  est. frequency of background counts (Hz)
- gvx  Galvo center position in x
- gvy  .. in y
- itr  iteration id
- lcx  localization position rel. to beam (m) in x
- lcy  .. in y
- lcz  .. in z
- lnc  loc non corrected 
- loc  absolute localization precision in x,y,z (m), NaN means
  no localization possible/done. Localization corrected with beam line correaction (for newer system)
  If lnc =/= loc then data has been corrected. 
- sta  debug parameter (state, number that encodes the abort
  criteria)
- tic  timestamp of every iteration / FPGA ticks 40MHz = (25ns)

All fields are NxM vectors for M iterations. Some are NxMx3: ext,
loc, lnc mostly for x,y,z.

Notes:

- fbg is not very often updated, does not contain many unique values,
  accuracy is doubtful
- background detection can be turn on/off with "BgcSense" true or
  false in sequence)
- eco, ecc, efo, efc may be background corrected (subtracted fbg)
  but there are hints that they might not be fully corrected
  (fbg determined wrongly)

### Number of photons
R.S. email:
`` For the first localization of a trace, you need the combined number of the “phtLimit” parameters (of all iterations of your sequence) as a lower limit for total detected photons. After the first localization succeeded, additional localizations require only as much photons as needed for iterations from “headstart” on to the final one.``

the parameter 'phtLimit' and 'headstart' can be found in the zarr directory 'mfx/.zattrs'. 
The 'headstart' is negative? The 'phtLimit' appears several time in the metadata with different values. 

### Beads for beam adjustment

You access the bead positions via the Zarr archive in the data directory (usually C:\Data\<data uid>\zarr) 
while the measurement is open in Imspector. Within the archive structure, their localization data 
is stored at /grd/mbm/R<nn>.

### Imspector/Parview 
Save data of paraview to file to load different measurements in one sample.
From Imspector launch paraview and than save data. 

# Loadig of msr files using python

Impsector comes with python libraries for the handling of the files. 
specpy is located at `C:\Imspector\Versions\IDOFIMSPECTORVERSION\python\PythonVersion-NumpyVersion`
e.g. `C:\Imspector\Versions\16.3.15620-m2205-win64-MINFLUX_BASE\python\Python3.9.2-NumPy1.20.1`
create environment with appropriate python version

```
conda create -n imspector-env python=3.9.2
conda activate imspector
cd C:\Imspector\Versions\16.3.15620-m2205-win64-MINFLUX_BASE\python\Python3.9.2-NumPy1.20.1
pip install specpy-1.2.3-cp39-cp39-win_amd64.whl
```


## load a msr 
`file=specpy.File("name.msr")` 
list the contained Minflux data sets (see dir() for the respective functions, they should be recognizable by name). 
Then, use unpack() to extract the data as a zarr to disk. 
This way, you can access all the mbm reference bead infos under /zarr/grd/mbm as arrays. 
The names of the beads that were used for the lnc->loc conversion are stored as "used" attribute of the mbm group.



# Python install instructions
Install conda or miniconda (https://docs.conda.io/en/latest/miniconda.html). Start anaconda prompt and install in the base environment

```
conda install -c conda-forge jupyter
conda install -c conda-forge nb_conda_kernels
``` 

The second package is to get see the different conda environments in jupyter automatically. Otherwise you need to do 
  
```
conda activate myenv
python -m ipykernel install --user --name myenv --display-name "Python (myenv)"
```
  

Create a working environment, we name it `mfx`

```
conda create -n mfx python=3.9 
conda activate mfx
```

## Use the yml file
Use the provided yml file and install all the needed packages
  
```
conda env create -f mfx.yml
``
  
Install imspector specific package

```
cd C:\Imspector\Versions\16.3.15620-m2205-win64-MINFLUX_BASE\python\Python3.9.2-NumPy1.20.1
pip install specpy-1.2.3-cp39-cp39-win_amd64.whl
```

Install the current python code from this repo. Go to the directory of this repo `python` type

```
pip install -e mfx
```
This is an editable install (the -e ). Any new update of the package will be taken into account and be accessible. 
  
## python packages in detail
  
Tools for the data processing

```
conda install -c conda-forge scikit-learn 
conda install -c conda-forge matplotlib
conda install -c conda-forge pandas
```


Zarr data handling
```
conda install -c conda-forge zarr
```

VTK stuff
```
conda install -c conda-forge vtk
conda install -c conda-forge pyevtk
```

Jupyter embedding
```
conda install ipykernel
```

Imspector specific package

```
cd C:\Imspector\Versions\16.3.15620-m2205-win64-MINFLUX_BASE\python\Python3.9.2-NumPy1.20.1
pip install specpy-1.2.3-cp39-cp39-win_amd64.whl
```

Finally install the python package in an editable version (still alpha). Go to the directory of this repo `python` type

```
pip install -e mfx
```

Any new update of the package will be taken into account and be accessible. 

# Start the software
Start a anaconda prompt. Make sure you are in the base environment. Eventually do ``` conda deactivate``` to get out of another environment. 
Go to the directory containing the code. Start a anacondo prompt and type

```
jupyter notebook
```


