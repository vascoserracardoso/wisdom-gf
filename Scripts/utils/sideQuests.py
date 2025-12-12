from astropy.io import fits
import os
from ..analysis import Plotter

# --- Get run specific values class ---
class getValues:
    def __init__(self):
        pass

    def lowestResolution(galaxies, max_keys):
        resolutions_list = []
        for galaxy_name, subdict in galaxies.items():
            if len(subdict) == max_keys:
                D = galaxies[galaxy_name]["distance"] 

                file_name = galaxy_name + '_mom0.fits'
                infile = './Galaxies/' + galaxy_name + '/FITS/inputs/' + file_name
                with fits.open(infile) as hdul:
                    hdr  = hdul[0].header.copy()
                ddeg1 = abs(hdr['CDELT1'])
                dtheta1 = ddeg1 * 3600.0 
                nx = hdr['NAXIS1']
                ny = hdr['NAXIS2']
                resolution = D * dtheta1 / 206265.0
                resolutions_list.append([galaxy_name, resolution, nx, ny]) # saved galaxy name and resolution
        lowest_resolution = max(resolutions_list, key=lambda x: x[1]) # lambda is an anonymous function returning the second element (resolution)

        return lowest_resolution
    
    def get_folder_size(path):
        total_size = 0
        for root, dirs, files in os.walk(path):
            for file in files:
                fp = os.path.join(root, file)
                # Skip broken symlinks
                if os.path.isfile(fp):
                    total_size += os.path.getsize(fp)
        size_mbytes = total_size / 1048576 # Convert to megabytes
        print(f"Total size of {path}: {size_mbytes:.1f} MB")

    def manual_center(self, gas, star, show_figures):
        center = input("Enter center coordinates as 'x,y': ")
        if center:
            center_x, center_y = center.split(',')
            center_x, center_y = float(center_x), float(center_y)
            plotter = Plotter(show_figures)
            plotter.plot_center(gas, center_x, center_y)
            plotter.plot_center(star, center_x, center_y)
            again = input("Again? (y/n): ")
            if again.lower() == 'y':
                self.manual_center(gas, star, show_figures)

            return center_x, center_y
        
        else: 
            return None