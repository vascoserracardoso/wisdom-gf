#!/usr/bin/env python3

# --- Import packages ---
import numpy as np
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils.segmentation import detect_sources
from astropy.wcs import WCS
from scipy.ndimage import zoom
import warnings
from astropy.wcs.wcs import FITSFixedWarning
from astropy.convolution import Gaussian2DKernel, convolve_fft

# --- Import custom modules ---
from .Dictionary import galaxies
from .analysis import Plotter
from .utils import userChoice

# --- Command line interface ---
def cli():
    imported_mom0 = False
    options = galaxies.keys()
    galaxy_name = userChoice(options).galaxyChoice()
    show_figures = input("Show figures? (y/n): ").lower() == 'y'
    main(galaxy_name,
         imported_mom0,
         show_figures,
         match_maps=False,
         lowest_resolution=None)
    

def main(galaxy, imported_script, show_figures, match_maps, lowest_resolution): 

    imported = imported_script

    # --- Import moment-0 map  and load data & header ---
    print("Importing moment-0 map...")
    file_name = galaxy + '_mom0.fits'
    infile = './Galaxies/' + galaxy + '/FITS/inputs/' + file_name
    print("- Input file:", infile)
    with fits.open(infile) as hdul:
        hdr  = hdul[0].header.copy()
        mom0 = hdul[0].data.astype(float)
    print(f"- mom0 grid computed (shape: {mom0.shape})")
    try:
        print(f"  Quantity (units): {hdr['BTYPE']} ({hdr['BUNIT']})")
    except KeyError as e:
        print(f"  <Missing header keyword: {e}>")
        print("  Assuming Integrated Intensity in units Jy/beam km/s")
    ddeg1 = abs(hdr['CDELT1']) # degrees
    ddeg2 = abs(hdr['CDELT2'])
    dtheta1 = ddeg1 * 3600.0   # arcsec
    dtheta2 = ddeg2 * 3600.0 
    print(f"- Pixel size: {dtheta1:.3f} arcsec in RA, {dtheta2:.3f} arcsec in DEC")
    D = galaxies[galaxy]["distance"]
    print(f'- Using distance to the galaxy D = {D:.2e} pc')
    Dtheta1 = D * dtheta1 / 206265.0 # pc
    Dtheta2 = D * dtheta2 / 206265.0
    print(f"- Pixel size: {Dtheta1:.3e} pc in RA, {Dtheta2:.3e} pc in DEC")

    # --- Inicialize plotter with the value of show_figures ---
    plotter = Plotter(show_figures)
    plotter.fits_image(mom0, file_name)

    # --- Create mask, segmentation and identification of main segment ---
    print("\nCreating mask...")
    momclip = hdr['MOMCLIP'] # 3σ noise in Jy/beam km/s
    if momclip is None:
        raise KeyError("MOMCLIP keyword not found in header")
    print(f"- 3\u03C3 noise level (MOMCLIP): {momclip:.3e} Jy/beam km/s")
    segm = detect_sources(mom0, threshold=momclip, npixels=10)
    momclip_2 = momclip
    main_id = np.argmax(segm.areas)
    mask = segm.data == (main_id + 1)
    mom0_masked_segmented = mom0 * mask
    print("- Identified main segment")
    plotter.fits_image(mask, file_name + " masked")

    # --- Calculate H2 mass surface density Σ_H2 ---
    print("\nCalculating Σ_H2...")
    # Determine line transition from rest frequency
    nu_rest = hdr['RESTFRQ'] * u.Hz
    if (nu_rest < 116.0e9 * u.Hz) and (nu_rest > 115.0e9 * u.Hz):
        print(f"- Rest frequency of CO(1-0): {nu_rest:.3e}")
        print("  Integrated intensity comes as I_10")
        mom0_converted = mom0
        mom0_masked_segmented_converted = mom0_masked_segmented
        momclip_2 = momclip_2
        print("  No need to use line ratio conversion factor")
    elif (nu_rest < 231.0e9 * u.Hz) and (nu_rest > 230.0e9 * u.Hz):
        print(f"- Rest frequency of CO(2-1): {nu_rest:.3e}")
        print("  Integrated intensity comes as I_21")
        R21 = 0.65
        mom0_converted = mom0 / R21
        mom0_masked_segmented_converted = mom0_masked_segmented / R21
        momclip_2 = momclip_2 / R21
        print(f"  I_10 calculated from line ratio conversion factor R21: {R21:.3f} (Leroy et al. 2021)")
    elif (nu_rest < 346.0e9 * u.Hz) and (nu_rest > 345.0e9 * u.Hz):
        print(f"- Rest frequency of CO(3-2): {nu_rest:.3e}")
        print("  Integrated intensity comes as I_32")
        R31 = 0.31 # CO(3-2)/CO(1-0) ratio, typical value for ETGs (Davis et al. 2014)
        mom0_converted = mom0 / R31
        mom0_masked_segmented_converted = mom0_masked_segmented / R31
        momclip_2 = momclip_2 / R31
        print(f"  I_32 calculated from line ratio conversion factor R31: {R31:.3f} (Leroy et al. 2021)")
    rad1 = ddeg1 * u.deg.to(u.rad)  # degrees to radians
    rad2 = ddeg2 * u.deg.to(u.rad)
    cellsize = (rad1 * rad2) # pixel area in steradians
    print(f"- Pixel area in steradians: {cellsize:.3e}")
    # Beam solid angle (steradians), BMAJ, BMIN are in degrees in the header
    bmaj = hdr['BMAJ'] * u.deg
    bmin = hdr['BMIN'] * u.deg
    omega_beam = (np.pi * bmaj.to(u.rad) * bmin.to(u.rad) / (4.0 * np.log(2.0)))  # → steradian
    print(f"- Beam solid angle (Ω_beam): {omega_beam:.3e} (steradians)")
    # Multiply by 1 Jy = W m⁻² Hz⁻¹ = J s⁻¹ m⁻² Hz⁻¹; divide by Ω_beam → K
    factor_K_per_Jy = u.Jy.to(u.K, equivalencies=u.brightness_temperature(nu_rest, beam_area=omega_beam))
    print(f"- Unit Conversion factor: {factor_K_per_Jy:.3e} K per (Jy/beam)")
    # Apply to moment-0 map (Jy/beam km/s → K km/s)
    Ico = mom0_converted * factor_K_per_Jy
    Ico_masked = mom0_masked_segmented_converted * factor_K_per_Jy
    momclip_2 = momclip_2 * factor_K_per_Jy
    print(f"- I_CO map computed in K km/s")
    alphaCO = galaxies[galaxy]["alpha_CO"] # 4.3 is the Milky Way value
    print(f"- Using α_CO = {alphaCO} M_sun/(K km/s pc^2)")
    Sigma = Ico * alphaCO * (u.Msun / u.pc**2) # M_sun/pc^2
    Sigma_masked = Ico_masked * alphaCO 
    momclip_2 = momclip_2 * alphaCO
    print("- Σ_H2 map computed in M_sun/pc^2")

    plotter.fits_image(Sigma_masked, galaxy + " - Σ_H2 masked")

    Sigma_masked = Sigma_masked * (u.Msun / u.pc**2) # M_sun/pc^2

    # --- Update headers and write outputs ---
    outfileSigmaH2 = './Galaxies/' + galaxy + '/FITS/outputs/' + galaxy + '_Sigma_H2.fits'
    hdr['BUNIT'] = 'M_sun/pc^2'
    hdr.add_history(f"Sigma[H2] = {Sigma.mean():.3e} M_sun/pc^2 from N[H2]")
    fits.writeto(outfileSigmaH2, Sigma.value, hdr, overwrite=True)
    print(f"- Wrote Σ_H2 map in M_sun/pc^2 to: {outfileSigmaH2}")
    # Masked output
    outfileSigmaH2Mask = './Galaxies/' + galaxy + '/FITS/outputs/' + galaxy + '_Sigma_H2_mask.fits'
    hdr['BUNIT'] = 'M_sun/pc^2'
    hdr.add_history(f"Sigma[H2] = {Sigma_masked.mean():.3e} M_sun/pc^2 from N[H2]")
    fits.writeto(outfileSigmaH2Mask, Sigma_masked.value, hdr, overwrite=True)
    print(f"- Wrote Σ_H2 masked map in M_sun/pc^2 to: {outfileSigmaH2Mask}")

    # Convert angular -> physical (pc/pix): pc = D * (arcsec / 206265)
    arcsec_per_deg = 3600.0
    pc_pix_x = (D * (ddeg1 * arcsec_per_deg) / 206265.0)
    pc_pix_y = (D * (ddeg2 * arcsec_per_deg) / 206265.0)

    # --- Match all maps to the lowest physical resolution (pc) - Beam matching and pixel matching ---
    if (match_maps and lowest_resolution[0] != galaxy) or (match_maps and lowest_resolution[1] != pc_pix_x):
        # Initialize reference file
        gas_file_ref  = './Galaxies/' + lowest_resolution[0] + '/FITS/outputs/' + lowest_resolution[0] + '_Sigma_H2_mask.fits'     # Molecular (H2) gas surface-density map
        # Load & squeeze to 2D
        with fits.open(gas_file_ref) as hdul:
            gas_hdu_ref  = hdul[0]
            gas_data_ref = gas_hdu_ref.data.squeeze().astype(float)
            gas_hdr_ref  = gas_hdu_ref.header
            with warnings.catch_warnings(): # Ignore FITSFixedWarning
                warnings.simplefilter('ignore', FITSFixedWarning) # Ignore WCS warnings
                gas_wcs_ref = WCS(gas_hdr_ref).celestial
        # Beam matching
        if 'BMAJ' in gas_hdr_ref and 'BMIN' in gas_hdr_ref:
            gas_bmaj  = gas_hdr_ref['BMAJ'] * 3600.0    # arcsec
            pix_scale = abs(gas_hdr_ref['CDELT1']) * 3600.0
            beam_pc = (hdr['BMAJ'] * 3600 * galaxies[galaxy]['distance'] / 206265) # only bmaj for the following calculation since the beams are usually circularized
            beam_ref_pc = (gas_hdr_ref['BMAJ'] * 3600 * galaxies[lowest_resolution[0]]['distance'] / 206265)
            # It only makes sense to use convolution in beams that have lower angular resolution
            if (beam_pc < beam_ref_pc) and (hdr['BMAJ'] * 3600 < gas_bmaj): 
                fwhm_diff = np.sqrt(gas_bmaj**2 - (hdr['BMAJ']*3600)**2)
                sigma_pix = fwhm_diff / (2.355 * pix_scale)
                kernel    = Gaussian2DKernel(sigma_pix)
                Sigma_mask_conv = convolve_fft(Sigma_masked, kernel, normalize_kernel=True)
                print(f"- Map convolved to {lowest_resolution[0]}'s beam.")
                K = kernel.array / kernel.array.sum() # propagation of error over convolution
                conv_factor = np.sqrt((K**2).sum())
                momclip_2 = momclip_2 * conv_factor
            else:
                print("- Beam already coarser than reference (physically or angularly); skipping convolution.")
                Sigma_mask_conv = Sigma_masked
        else:
            print("No gas beam info available; skipping convolution.")
            # no gas beam info → leave Sigma_H2 untouched
            Sigma_mask_conv = Sigma_masked

        # --- Pixel matching - Resampling to target physical pixel size and shape ---
        print(f"\nResampling Σ_H2 masked to physical pixel size of {lowest_resolution[1]} pc/pix ({lowest_resolution[0]}) and shape {lowest_resolution[2]}x{lowest_resolution[3]}...")
        # Target physical pixel size (pc/pix) and target shape (nx, ny)
        target_pc_per_pix = float(lowest_resolution[1])  # pc/pix (scalar)
        nx_target = int(lowest_resolution[2])
        ny_target = int(lowest_resolution[3])
        ref_shape = (ny_target, nx_target)
        # Convert angular -> physical (pc/pix): pc = D * (arcsec / 206265)
        arcsec_per_deg = 3600.0
        pc_per_pix_x = (D * (ddeg1 * arcsec_per_deg) / 206265.0)
        pc_per_pix_y = (D * (ddeg2 * arcsec_per_deg) / 206265.0)
        print(f"- Current pc/pix (x,y): {pc_per_pix_x:.4f}, {pc_per_pix_y:.4f}")
        print(f"- Target pc/pix: {target_pc_per_pix:.4f}")
        # Zoom factors for scipy.zoom: new_size = old_size * zoom_factor
        zoom_x = pc_per_pix_x / target_pc_per_pix
        zoom_y = pc_per_pix_y / target_pc_per_pix
        print(f"- Zoom factors (y, x): {zoom_y:.4f}, {zoom_x:.4f}  (new_n = old_n * zoom)")
        # Resample (order=1 bilinear; choose order=0/3 depending on desired accuracy/performance)
        resampled = zoom(Sigma_mask_conv.value, (zoom_y, zoom_x), order=1, prefilter=False)
        print(f"- Resampled shape: {resampled.shape}")
        # Center resampled into target array (pad with nan or crop symmetrically)
        target = np.full(ref_shape, np.nan, dtype=resampled.dtype)
        ny_t, nx_t = ref_shape
        ny_r, nx_r = resampled.shape
        start_y = (ny_t - ny_r) // 2
        start_x = (nx_t - nx_r) // 2
        if start_y >= 0 and start_x >= 0:
            # resampled smaller -> pad
            target[start_y:start_y + ny_r, start_x:start_x + nx_r] = resampled
        else:
            # resampled larger -> crop symmetrically
            crop_y0 = max(0, -start_y)
            crop_x0 = max(0, -start_x)
            crop_y1 = crop_y0 + min(ny_t, ny_r)
            crop_x1 = crop_x0 + min(nx_t, nx_r)
            target[:] = np.nan
            ty0 = max(0, start_y)
            tx0 = max(0, start_x)
            ty1 = ty0 + (crop_y1 - crop_y0)
            tx1 = tx0 + (crop_x1 - crop_x0)
            target[ty0:ty1, tx0:tx1] = resampled[crop_y0:crop_y1, crop_x0:crop_x1]
        # Replace NaN or inf with 0
        target = np.nan_to_num(target, nan=0.0, posinf=0.0, neginf=0.0)
        print(f"- Final target shape: {target.shape}  finite pixels: {np.isfinite(target).sum()}")
        # theta_rad = target_pc_per_pix / D  -> radians; convert to degrees
        theta_rad = (target_pc_per_pix / D)
        deg_per_pix_ref = np.degrees(theta_rad)
        # Preserve original sign/orientation for axis 1 if present
        sign_x = np.sign(hdr.get('CDELT1', 1.0))
        sign_y = np.sign(hdr.get('CDELT2', 1.0))

        # --- Update header: set CDELT to match physical target (convert pc/pix -> deg/pix) ---
        new_hdr = hdr.copy()
        new_hdr['CDELT1'] = sign_x * deg_per_pix_ref
        new_hdr['CDELT2'] = sign_y * deg_per_pix_ref
        # Set CRPIX to the center of the target array (pixel coordinates are 1-based in FITS headers)
        new_hdr['CRPIX1'] = nx_target / 2.0
        new_hdr['CRPIX2'] = ny_target / 2.0
        # Optionally update NAXIS1/2 to the new shape
        new_hdr['NAXIS1'] = nx_target
        new_hdr['NAXIS2'] = ny_target
        new_hdr['BMAJ'] = gas_hdr_ref['BMAJ']
        new_hdr['BMIN'] = gas_hdr_ref['BMIN']
        if 'BPA' in gas_hdr_ref:
            new_hdr['BPA'] = gas_hdr_ref['BPA']

        plotter.fits_image(target, galaxy + " - Σ_H2 masked resampled (phys match)")

        # --- Mask to clear out points created due to the blurring - Using the noise processed up to convolution, because during regreding there's no significant chance
        segm = detect_sources(target, threshold=momclip_2, npixels=1) 
        # Identify largest segment
        main_id = np.argmax(segm.areas)
        mask = segm.data == (main_id + 1)
        target_masked = target * mask
        print("- Identified main segment")
        plotter.fits_image(target_masked, galaxy + " - Σ_H2 masked resampled (phys match) masked")

        # --- Save result ---
        outfile = f'./Galaxies/{galaxy}/FITS/outputs/{galaxy}_Sigma_H2_mask_resampled_to_{lowest_resolution[0]}.fits'
        fits.writeto(outfile, target_masked, new_hdr, overwrite=True)
        print(f"- Saved resampled (physical-match) Σ_H2 masked (twice) map to: {outfile}")

    # --- Calculate Total Mass ---
    print("\nCalculating total mass in the map...")
    A_pix = Dtheta1 * Dtheta2  # pc^2
    print(f"- Pixel area: {A_pix:.3e} pc²")
    M_H2 = (Sigma.sum()) * A_pix
    print(f"- Total mass in the map: {M_H2:.3e}")
    # Mass in the masked region
    M_H2 = Sigma_masked.sum() * A_pix * u.pc**2 # So it spits out solMass
    print(f"- Total mass in masked region: {M_H2:.3e}")

    # --- Write M_H2 to PhysicalQuantitiesOutputs.txt ---
    if imported:
        with open(f'./Galaxies/{galaxy}/PhysicalQuantitiesOutputs.txt', 'a') as f:
            f.write(f"M_H2                {M_H2:.6e}\n")
        print(f"- Saved M_H2 to ./Galaxies/{galaxy}/PhysicalQuantitiesOutputs.txt")


if __name__ == "__main__":
    cli()


