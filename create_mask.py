import numpy as np
from netCDF4 import Dataset

# Open the input NetCDF file
input_file = '../Moorings_2007m01.nc'  # Change this to your input file
output_file = 'FullArctic.nc'  # Change this to your desired output file

# Open the dataset
with Dataset(input_file, 'r') as ds:
    # Read the variable (change 'sic' to your variable name)
    sic = ds.variables['sic'][:][0,:,:]
    # Get latitude and longitude (assuming they are named 'latitude' and 'longitude')
    lat = ds.variables['latitude'][:]
    lon = ds.variables['longitude'][:]
   
    print(np.shape(sic))
    # Create the mask
    mask = np.zeros_like(sic, dtype=np.float64)  # Initialize the mask with zeros
    mask = np.zeros(np.shape(sic))
    print(mask[0,0])

    # Extract _FillValue attribute
    fill_value = ds.variables['sic']._FillValue if '_FillValue' in ds.variables['sic'].ncattrs() else None

    # Set mask to 1 where sic is between 0 and 1
    mask[(sic >= 0) & (sic <= 1)] = 1
    print(mask[0,0])
    
    # Handle _FillValue correctly
    if fill_value is not None:
        mask[sic == fill_value] = 0  # Set mask to 0 where sic equals _FillValue
    print(mask[0,0])

    # Write the output to a new NetCDF file
    with Dataset(output_file, 'w', format='NETCDF4') as ds_out:
        # Create dimensions
        ds_out.createDimension('y', ds.dimensions['y'].size)        # Y dimension
        ds_out.createDimension('x', ds.dimensions['x'].size)        # X dimension
    
        # Create variables in the new dataset
        lat_var = ds_out.createVariable('latitude', lat.dtype, ('y','x')) 
        lon_var = ds_out.createVariable('longitude', lon.dtype, ('y','x',))  
        mask_var = ds_out.createVariable('mask', mask.dtype, ('y', 'x'))  
    
        # Assign data to variables
        lat_var[:,:] = lat
        lon_var[:,:] = lon
        mask_var[:, :] = mask  # Note: Adjust the dimensions accordingly if needed
    
        # Copy attributes from the original dataset (if needed)
        # mask_var.setncatts(ds.variables['sic'].__dict__)  # Copy attributes from original variable

print(f'Mask created and saved to {output_file}')
