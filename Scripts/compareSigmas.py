#!/usr/bin/env python3

# --- Import packages ---
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from reproject import reproject_interp
from astropy.convolution import Gaussian2DKernel, convolve_fft
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterExponent
from matplotlib.colors import Normalize
import warnings
from astropy.wcs.wcs import FITSFixedWarning

# --- Import custom modules ---
from .utils import getValues
from .analysis import Plotter, set_rc_params
from .Dictionary import galaxies
from .utils import userChoice


# --- Command line interface ---
def cli():
    imported_compare = False
    options = galaxies.keys()
    galaxy_name = userChoice(options).galaxyChoice()
    show_figures = input("Show figures? (y/n): ").lower() == 'y'
    main(galaxy_name,
         imported_compare,
         show_figures,
         match_maps=False,
         lowest_resolution=None)

def main(galaxy, imported_script, show_figures, match_maps, lowest_resolution):
    set_rc_params(0.5) # setting MNRAS style with a multiplier of 0.8 for font sizes

    imported = imported_script

    # --- Import file and load header ---
    print("Importing input files and regridding...")
    file = './Galaxies/' + galaxy + '/FITS/inputs/' + galaxy + '_mom0.fits'
    file_hdr = fits.getheader(file)
    D = galaxies[galaxy]["distance"]
    # Current angular pixel sizes (deg/pix)
    cdelt_x_deg = abs(file_hdr['CDELT1'])  # deg/pix along X (RA)
    cdelt_y_deg = abs(file_hdr['CDELT2'])  # deg/pix along Y (Dec)
    # Convert angular -> physical (pc/pix): pc = D * (arcsec / 206265)
    arcsec_per_deg = 3600.0
    pc_pix_x = (D * (cdelt_x_deg * arcsec_per_deg) / 206265.0)
    pc_pix_y = (D * (cdelt_y_deg * arcsec_per_deg) / 206265.0)

    # --- Match maps option ---
    if (match_maps and lowest_resolution[0] != galaxy) or (match_maps and lowest_resolution[1] != pc_pix_x):
        gas_file      = './Galaxies/' + galaxy + '/FITS/outputs/' + galaxy + '_Sigma_H2_mask_resampled_to_' + lowest_resolution[0] + '.fits'     # Molecular (H2) gas surface-density map
        star_file     = './Galaxies/' + galaxy + '/FITS/outputs/' + galaxy + '_Sigma_star_resampled_to_' + lowest_resolution[0] + '.fits'   # Stellar surface-density map
    else:
        gas_file      = './Galaxies/' + galaxy + '/FITS/outputs/' + galaxy + '_Sigma_H2_mask.fits'     # Molecular (H2) gas surface-density map
        star_file     = './Galaxies/' + galaxy + '/FITS/outputs/' + galaxy + '_Sigma_star.fits'   # Stellar surface-density map
    
    # --- Load & squeeze to 2D ---
    with fits.open(gas_file) as hdul:
        gas_hdu  = hdul[0]
        gas_data = gas_hdu.data.squeeze().astype(float)
        gas_hdr  = gas_hdu.header
        with warnings.catch_warnings(): # Ignore FITSFixedWarning
            warnings.simplefilter('ignore', FITSFixedWarning) # Ignore WCS warnings
            gas_wcs = WCS(gas_hdr).celestial
    with fits.open(star_file) as hdul:
        star_hdu  = hdul[0]
        star_data = star_hdu.data.squeeze().astype(float)
        star_hdr  = star_hdu.header
        with warnings.catch_warnings(): # Ignore FITSFixedWarning
            warnings.simplefilter('ignore', FITSFixedWarning) # Ignore WCS warnings
            star_wcs  = WCS(star_hdr).celestial
    #gas_wcs = WCS(gas_hdr, naxis=2)
    #star_wcs = WCS(star_hdr, naxis=2)
    ny, nx = gas_data.shape
    #assert (ny, nx) == star_data.shape, "Gas and star data must have same shape!"
    print("Gas data shape:", gas_data.shape)
    print("Stellar data shape:", star_data.shape)


    # --- In case the user does not want to match maps ---
    if match_maps is not True:
        # --- Reproject stellar onto gas grid ---
        star_reproj, footprint = reproject_interp(
            (star_data, star_wcs),
            gas_wcs,
            shape_out=(ny, nx)
        )
        print("Gas data shape:", star_reproj.shape)
        print("Reprojected stellar shape:", star_reproj.shape)
        print("- Stellar map onto gas grid using header information")

        # --- Beam matching ---
        print("Maximum stellar intensity before beam matching:", star_reproj.max())
        # Beam matching: always convolve star -> gas beam ---
        if 'BMAJ' in gas_hdr and 'BMIN' in gas_hdr:
            gas_bmaj  = gas_hdr['BMAJ'] * 3600.0    # arcsec
            pix_scale = abs(gas_hdr['CDELT1']) * 3600.0
            # star beam = 0, so fwhm_diff = gas_bmaj
            fwhm_diff = gas_bmaj
            sigma_pix = fwhm_diff / (2.355 * pix_scale)
            kernel    = Gaussian2DKernel(sigma_pix)
            star_conv = convolve_fft(
                star_reproj, kernel,
                normalize_kernel=True
                )
        else:
            print("No gas beam info available; skipping convolution.")
            # no gas beam info → leave star untouched
            star_conv = star_reproj
        print("- Stellar map convolved to gas beam")
        print("Maximum stellar intensity after beam matching:", star_conv.max())

    gal_center_x = nx//2
    gal_center_y = ny//2
    print(f"- Using geometrical center of the image as center of the object: ({gal_center_x}, {gal_center_y})")
    if match_maps is True:
        star_conv = star_data # In case of no beam matching, use original star data
    plotter = Plotter(show_figures)
    plotter.plot_center(gas_data,gal_center_x,gal_center_y)
    plotter.plot_center(star_conv,gal_center_x,gal_center_y)

    # --- Deproject to face-on ---
    inclination   = galaxies[galaxy]["inclination"] # deg;
    print("- Using inclination:", inclination), " to turn maps into Face-on"
    inc_rad = np.radians(inclination)
    gas_data  = gas_data  / np.cos(inc_rad) #Face-on gas map
    star_conv = star_conv  / np.cos(inc_rad)

    # --- Save outputs ---
    if match_maps is not True:
        print("\nSaving outputs...")
        rep_beam_file = './Galaxies/' + galaxy + '/FITS/outputs/' + galaxy + '_Sigma_star_rep_beam.fits'
        hdr_star = gas_hdr.copy()
        hdr_star['BUNIT']   = star_hdr.get('BUNIT', 'M_sun/pc^2')
        hdr_star['HISTORY'] = 'Reprojected & beam-matched stellar map'
        fits.writeto(rep_beam_file, star_conv, hdr_star, overwrite=True)
        print(f"- Wrote reprojected & beam-matched stellar map to: {rep_beam_file} (done according to its gas map correspondent)")

    # --- Calculate map of distances to the center of the galaxy ---
    print("\nCalculating map of distances to the center of the galaxy...")
    yy, xx = np.indices(gas_data.shape)
    # Distance from galaxy center (in pixels)
    distances = np.sqrt((xx - gal_center_x)**2 + (yy - gal_center_y)**2)


# --- Plot pixel-by-pixel comparison ---
    print("\nPlotting pixel-by-pixel comparison...")
    # Parameters for the best fit in Lin et al. (2020)
    a = 1.10
    b = -1.95
    sig = 0.20 # Probably the uncertainty of a
    ro = 0.76 # Pearson correlation coefficient
    print(f'- Using Lin et al. (2020) fit: a = {a:.2f}, b = {b:.2f}, σ = {sig:.2f}, ρ = {ro:.2f}')
    # Filter out non-positive values
    mask  = (gas_data > 0) & (star_conv > 0)
    print(f"- Masked {np.sum(~mask)} pixels with non-positive values")
    # Flatten data
    x = star_conv[mask].ravel() * 10**6 # from pc² to kpc²
    y = gas_data[mask].ravel() * 10**6 # from pc² to kpc²
    # Stack x and y as columns
    xy = np.column_stack((x, y))
    # Save to file with two columns: x and y
    np.savetxt(f'./Galaxies/{galaxy}/rMGMSOutputs.txt', xy, header='Sigma_star[Sol_mass kpc^-2] Sigma_H2[Sol_mass kpc^-2]', fmt='%.6e')
    print(f"- Saved Sigma_H2 vs Sigma_star to ./Galaxies/{galaxy}/rMGMSOutputs.txt")
    print("- Flattened data for plotting (length: {})".format(len(x)))
    # Generate values according to x and y
    x_lin = np.logspace(np.log10(1e5), np.log10(1e12), 100)
    y_lin = 10**b * x_lin**a
    sig_mult = 1
    y_lin_upper = y_lin * 10 ** (sig_mult*sig)
    y_lin_lower = y_lin * 10 ** (sig_mult*-sig)
    # Keep the same pixels that go into the scatter plot
    dist_pix = distances[mask].ravel() # 1-D, same length as x and y
    # --- optional physical units ---
    pixel_scale_arcsec = abs(gas_hdr['CDELT1']) * 3600 # deg → arcsec
    print(f"- Pixel scale: {pixel_scale_arcsec:.6f} arcsec/pixel")
    distance_mpc = galaxies[galaxy]["distance"] / 1e6  # Mpc
    kpc_per_arcsec = distance_mpc * 1e3 * np.radians(1/3600)   # small-angle approx
    print(f"- kpc per arcsec: {kpc_per_arcsec:.6f}")
    dist_for_colour = dist_pix * pixel_scale_arcsec * kpc_per_arcsec
    cbar_label = 'Projected radius (kpc)'
    norm = Normalize(vmin=dist_for_colour.min(), vmax=dist_for_colour.max())

    fig, ax = plt.subplots(figsize=(6,6))
    sc = ax.scatter(x, y,
                    c=dist_for_colour, cmap='plasma', norm=norm,
                    s=0.5, alpha=0.8, label=galaxy)
    # colour-bar
    fig.colorbar(sc, ax=ax, label=cbar_label)
    # set log scales
    ax.set_xscale('log')
    ax.set_yscale('log')
    # major ticks at powers of ten
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_locator(LogLocator(base=10))
    # formatter that only shows the exponent
    ax.xaxis.set_major_formatter(LogFormatterExponent(base=10))
    ax.yaxis.set_major_formatter(LogFormatterExponent(base=10))
    ax.set_xlabel(r'Log $\Sigma_{*}$ [M$_\odot$ kpc$^{-2}$]')
    ax.set_ylabel(r'Log $\Sigma_{\rm H_2}$ [M$_\odot$ kpc$^{-2}$]')
    ax.set_title('Pixel-by-Pixel: Σ_H2 vs. Σ_*')
    ax.grid(True, which='both', ls='--', alpha=0.3)
    # Plot Lin's fit
    ax.plot(x_lin, y_lin, color = 'black', lw = 1, label=f'Lin et al. (2020): a = {a}; b = {b}')
    ax.plot(x_lin, y_lin_upper, color = 'black', linestyle = '--', lw = 0.7, label = fr'$+\sigma$')
    ax.plot(x_lin, y_lin_lower, color = 'black', linestyle = '--', lw = 0.7, label = fr'$-\sigma$')
    ax.legend()
    fig.tight_layout()
    Lin_figure = f'./Galaxies/{galaxy}/Figures/' + galaxy + '_Lin2020_log.png'
    fig.savefig(Lin_figure, dpi=300)
    if show_figures:
        plt.show()
    plt.close()
    print("- Pixel-by-pixel comparison plot saved as:", Lin_figure)


if __name__ == "__main__":
    cli()


