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


### Beads for beam adjustment

You access the bead positions via the Zarr archive in the data directory (usually C:\Data\<data uid>\zarr) 
while the measurement is open in Imspector. Within the archive structure, their localization data 
is stored at /grd/mbm/R<nn>.

