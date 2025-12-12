#!/usr/bin/env python3

# --- Import packages ---
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from InquirerPy import inquirer
from scipy.ndimage import zoom

# --- Import custom modules ---
from .Dictionary import galaxies
from .utils import userChoice

# --- Command line interface ---
def cli():
    imported_mge = False
    options = galaxies.keys()
    galaxy_name = userChoice(options).galaxyChoice()
    show_figures = input("Show figures? (y/n): ").lower() == 'y'
    main(
        galaxy_name,
        imported_mge,
        show_figures,
        match_maps=False,
        lowest_resolution=None,
    )

def main(galaxy, imported_script, show_figures, match_maps, lowest_resolution):

    imported = imported_script

    # --- Calculate the stellar surface density map from MGE parameters ---
    print("Calculating Σ_* map from MGE parameters...")
    # Published MGE parameters
    Ij     = galaxies[galaxy]["Ij"]      # L_sol/pc^2
    sigmaj = galaxies[galaxy]["sigmaj"]  # arcsec
    qj     = galaxies[galaxy]["qj"]      # b/a
    print("- MGE parameters from North et al. (2019):")
    print(f"     Ij: {Ij}")
    print(f"     sigmaj: {sigmaj}")
    print(f"     qj: {qj}")
    # Load the gas map header
    sigmaH2_file = './Galaxies/' + galaxy + '/FITS/outputs/' + galaxy + '_Sigma_H2.fits'
    gas_hdr = fits.getheader(sigmaH2_file)

    # --- Build the Σ_* map on the gas map grid ---
    print("\nBuilding Σ_* map on gas map grid...")
    if match_maps is not True or lowest_resolution[0] == galaxy:
        deg_per_pix = abs(gas_hdr['CDELT1'])  # deg / pix
        pix_scale = deg_per_pix * 3600.0   # arcsec/pix
        ny = gas_hdr['NAXIS2']   # number of rows
        nx = gas_hdr['NAXIS1']   # number of columns
        x = (np.arange(nx) - nx//2) * pix_scale
        y = (np.arange(ny) - ny//2) * pix_scale
        X, Y = np.meshgrid(x, y)
        print("- Gas map header imported for grid construction")
        # Manually reconstruct the 2D Σ_L map, Σ_L(x,y) = Σ_j I_j * exp[-0.5 * (X^2 + (Y^2/q_j^2)) / σ_j^2]
        Sigma_L = np.zeros_like(X)
        for I, s, q in zip(Ij, sigmaj, qj):
            R2 = X**2 + (Y**2)/(q**2)
            Sigma_L += I * np.exp(-0.5 * R2 / s**2)
        # I is peak surface brightness
        print("- Σ_L map computed using the expression Σ_L(x,y) = Σ_j I_j * exp[-0.5 * (X^2 + (Y^2/q_j^2)) / σ_j^2]")
        print(f"  Quantity (units): Surface luminosity density (L_sun/pc^2)") # 1 solar luminosity is 3.828e26 W
        # Convert to stellar mass surface density via M/L
        ML = galaxies[galaxy]["ML"] # e.g. the paper’s best‐fit
        print("- Using M/L = {:.1f} conversion factor for stellar mass surface density".format(ML))
        Sigma_star = Sigma_L * ML  * u.M_sun / u.pc**2  # M_sun pc^-2
        print(f"- Σ_* grid computed (shape: {Sigma_star.shape})")

        # --- Construct a minimal 2D WCS ---
        ra_center = galaxies[galaxy]["ra_center"]
        dec_center = galaxies[galaxy]["dec_center"]
        w = WCS(naxis=2)
        w.wcs.crpix  = [gas_hdr['CRPIX1'], gas_hdr['CRPIX2']] # reference pixel
        #w.wcs.crpix  = [nx//2, ny//2]                         # center pixel
        w.wcs.crval  = [gas_hdr['CRVAL1'], gas_hdr['CRVAL2']] # RA, Dec centre (deg)
        w.wcs.cdelt = [gas_hdr['CDELT1'], gas_hdr['CDELT2']]  # match gas pixel scale
        w.wcs.ctype = [gas_hdr['CTYPE1'], gas_hdr['CTYPE2']]  # match projection
        w.wcs.cunit  = ["deg", "deg"]
        w.wcs.radesys = 'ICRS'
        w.wcs.equinox = 2000.0
        print(f"- Minimal WCS constructed with crpix = ({gas_hdr['CRPIX1']}, {gas_hdr['CRPIX2']}) and crval = ({gas_hdr['CRVAL1']}, {gas_hdr['CRVAL2']}) deg")

        # --- Update header and write output ---
        outfileSigmaStar = f'./Galaxies/{galaxy}/FITS/outputs/{galaxy}_Sigma_star.fits'
        hdr = w.to_header()
        hdr['BUNIT'] = 'M_sun/pc^2'
        hdr.add_history(f"Sigma_* from MGE (ML={ML}), grid 0.1\"/pix")
        fits.writeto(outfileSigmaStar, Sigma_star.value, hdr, overwrite=True)
        print(f"- Wrote Σ_* map in M_sun/pc^2 to: {outfileSigmaStar}")
        
        sigma_star_map = fits.getdata(outfileSigmaStar)
        plt.imshow(sigma_star_map, origin='lower', cmap='gray')

        # --- Plot objects with centers ---
        center_RA_DEC = SkyCoord(ra=ra_center*u.deg, dec=dec_center*u.deg, frame="icrs") # NED
        center_x, center_y = w.world_to_pixel(center_RA_DEC)
        print(f"- Using imported galaxy center: ({center_x}, {center_y}); ({ra_center}, {dec_center})")
        center_RA_DEC_gas = SkyCoord(ra=gas_hdr['CRVAL1']*u.deg, dec=gas_hdr['CRVAL2']*u.deg, frame="icrs") 
        center_x_gas, center_y_gas = w.world_to_pixel(center_RA_DEC_gas)
        print(f"- Using as reference pixel: ({center_x_gas}, {center_y_gas}); ({gas_hdr['CRVAL1']}, {gas_hdr['CRVAL2']})")
        plt.colorbar()
        plt.scatter(center_x, center_y, color='red', marker='.', s=25, label='Galaxy Center')
        plt.scatter(center_x_gas, center_y_gas, color='blue', marker='.', s=25, label='Gas Center')
        plt.title(galaxy + " - Σ_star.fits")
        plt.draw()
        if show_figures:
            plt.show()
        plt.close()


    D = galaxies[galaxy]["distance"]
    cdelt_x_deg = abs(gas_hdr['CDELT1'])  # deg/pix along X (RA)
    cdelt_y_deg = abs(gas_hdr['CDELT2'])  # deg/pix along Y (Dec)
    arcsec_per_deg = 3600.0
    pc_pix_x = (D * (cdelt_x_deg * arcsec_per_deg) / 206265.0)
    pc_pix_y = (D * (cdelt_y_deg * arcsec_per_deg) / 206265.0)

    # --- Regrid Σ_* to lowest resolution ---
    if (match_maps and lowest_resolution[0] != galaxy) or (match_maps and lowest_resolution[1] != pc_pix_x):
        print(f"\n--- Matching Σ_* map to reference grid ({lowest_resolution[0]})---")
        # Target resolution and size from lowest_resolution ---
        desired_resolution = float(lowest_resolution[1])  # pc / pixel
        nx = lowest_resolution[2]
        ny = lowest_resolution[3]
        dtheta_rad = (desired_resolution  / D) # radians / pix
        deg_per_pix = np.degrees(dtheta_rad)  # degrees / pix
        pix_scale = deg_per_pix * 3600.0  # arcsec / pix (for building X,Y)
        print(f"- Target grid: {nx} x {ny} pixels, {desired_resolution:.3f} pc/pix -> {deg_per_pix:.6e} deg/pix")
        # Build the new spatial grid in arcseconds centered on image center
        x = (np.arange(nx) - nx//2) * pix_scale
        y = (np.arange(ny) - ny//2) * pix_scale
        X, Y = np.meshgrid(x, y)
        print("- Constructed X,Y grid for new target (arcsec)")
        # Compute Sigma_L on new grid
        Sigma_L = np.zeros_like(X)
        for I, s, q in zip(Ij, sigmaj, qj):
            R2 = X**2 + (Y**2) / (q**2)
            Sigma_L += I * np.exp(-0.5 * R2 / s**2)
        print("- Σ_L computed on new grid")
        # Convert to Σ_*
        ML = galaxies[galaxy]["ML"]
        Sigma_star = Sigma_L * ML * u.M_sun / u.pc**2
        print(f"- Σ_* (shape {Sigma_star.shape}) computed")

        # --- Build a new minimal WCS with desired pixel scale & projection ---
        new_w = WCS(naxis=2)
        # Copy projection and coordinate types from gas header
        new_w.wcs.ctype = [gas_hdr['CTYPE1'], gas_hdr['CTYPE2']]
        new_w.wcs.cunit = ["deg", "deg"]
        new_w.wcs.radesys = getattr(gas_hdr, 'RADESYS', 'ICRS') if hasattr(gas_hdr, 'RADESYS') else 'ICRS'
        new_w.wcs.equinox = 2000.0
        # Use the gas map CRVAL (sky coordinate we want to preserve)
        crval_ra = gas_hdr['CRVAL1']
        crval_dec = gas_hdr['CRVAL2']
        new_w.wcs.crval = [crval_ra, crval_dec]
        # Set pixel scale (FITS convention: negative cdelt1 for RA axis if original uses that), else use [-deg_per_pix, deg_per_pix]
        cd1_sign = -1.0 if gas_hdr.get('CDELT1', -abs(deg_per_pix)) < 0 else 1.0
        new_w.wcs.cdelt = [cd1_sign * deg_per_pix, deg_per_pix]
        # Temporarily set crpix = 0 to find where the CRVAL would land; we'll then shift to desired pixel
        new_w.wcs.crpix = [0.0, 0.0]
        # The sky coordinate we want to anchor (use the gas header CRVAL)
        anchor_coord = SkyCoord(ra=crval_ra * u.deg, dec=crval_dec * u.deg, frame="icrs")
        # Where would that anchor land with current (crpix=0) new_w?
        px0, py0 = new_w.world_to_pixel(anchor_coord)  # float pixel positions
        # Desired pixel where that anchor should be located in the new image:
        # OPTION A (preserve original gas reference pixel location):
        desired_px = gas_hdr['CRPIX1']
        desired_py = gas_hdr['CRPIX2']
        # OPTION B (center the anchor in the new image instead):
        # desired_px = nx / 2.0
        # desired_py = ny / 2.0
        # Compute the required CRPIX so anchor_coord maps to desired pixel
        crpix_x = desired_px - px0
        crpix_y = desired_py - py0
        new_w.wcs.crpix = [crpix_x, crpix_y]
        print(f"- New WCS built. crval=({crval_ra:.6f}, {crval_dec:.6f}), crpix=({crpix_x:.3f}, {crpix_y:.3f}), cdelt=({new_w.wcs.cdelt[0]:.6e}, {new_w.wcs.cdelt[1]:.6e})")

        # --- Save to FITS ---
        outfileSigmaStar = f'./Galaxies/{galaxy}/FITS/outputs/{galaxy}_Sigma_star_resampled_to_' + lowest_resolution[0] + '.fits'
        hdr = new_w.to_header()
        # --- Update beam parameters (BMAJ, BMIN, BPA) ---
        # These stay in angular units, so no conversion to pc needed.
        # If the beam should remain the same physically, keep them identical.
        # Preserve the *same pixel-beam ratio*, scale them:
        scale_factor = abs(deg_per_pix / gas_hdr['CDELT1'])
        hdr['BMAJ'] = gas_hdr['BMAJ'] * scale_factor
        hdr['BMIN'] = gas_hdr['BMIN'] * scale_factor
        if 'BPA' in gas_hdr:
            hdr['BPA'] = gas_hdr['BPA']
        hdr['BUNIT'] = 'M_sun/pc^2'
        hdr.add_history(f"Sigma_* regridded to {desired_resolution:.3f} pc/pix ({nx}x{ny}); anchor CRVAL preserved")
        fits.writeto(outfileSigmaStar, Sigma_star.value, hdr, overwrite=True)
        print(f"- Wrote regridded Σ_* map to: {outfileSigmaStar}")

        # --- Plot the regridded Σ_* map ---
        sigma_star_map = fits.getdata(outfileSigmaStar)
        plt.imshow(sigma_star_map, origin='lower', cmap='gray')
        plt.colorbar()
        #plt.scatter(center_x, center_y, color='red', marker='.', s=25, label='Galaxy Center')
        #plt.scatter(center_x_gas, center_y_gas, color='blue', marker='.', s=25, label='Gas Center')
        plt.title(galaxy + " - Σ_star.fits")
        plt.draw()
        if show_figures:
            plt.show()
        plt.close()

    # --- Calculate total mass in the map ---
    print("\nCalculating total mass in the map...")
    D = D * u.pc
    A_pix = (D * pix_scale / 206265.0)**2
    print(f"- Pixel area in pc^2: {A_pix.to(u.pc**2):.3e} per pixel")
    M_star = (Sigma_star.sum()) * A_pix
    print(f"- Total mass in the map: {M_star:.3e}")

    # --- Write M_star to PhysicalQuantitiesOutputs.txt ---
    if imported:
        with open(f'./Galaxies/{galaxy}/PhysicalQuantitiesOutputs.txt', 'a') as f:
            f.write(f"M_star              {M_star:.6e}\n")
        print(f"- Saved M_star to ./Galaxies/{galaxy}/PhysicalQuantitiesOutputs.txt")

    # ------------------------------------------

if __name__ == "__main__":
    cli()



