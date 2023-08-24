== Prerequsites ==
NCIO: https://github.com/alex-robinson/ncio
datetime-fotran: https://github.com/wavebitscientific/datetime-fortran

== Data sets ==
ERA5 test: https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/fileServer/meomopendap/extract/SASIP/model-forcings/atmo_forcing/ERA5_Arctic/ERA5_msl_y2007.nc
Grid test: https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/fileServer/meomopendap/extract/SASIP/grids/NH_PS/25km_NH.nc

== RUNNING ON FRAM ==
The model is able to run on Fram with the following extra files:

- prerequisites from NCIO and datetime-fortran, as follows:
-- libncio.a, libdatetime.a, ncio_transpose.o, ncio_transpose.mod, ncio.o and ncio.mod from /cluster/projects/nn9878k/hregan/ABL/prerequisites/ncio
-- datetime_module.mod from /cluster/projects/nn9878k/hregan/ABL/prerequisites/datetime-fortran/build/include
These should be copied into the run directory

- a machine-specific make.macro.
-- there is a working one here: /cluster/projects/nn9878k/hregan/ABL/Earth_ABL_model/make.macro. In there, the tmp_dir should be set to something related to the user

- grid, initial and forcing files sufficient to run for 10 days are available here (all of the symbolic links):
/cluster/projects/nn9878k/hregan/ABL/Earth_ABL_model/data/ 


