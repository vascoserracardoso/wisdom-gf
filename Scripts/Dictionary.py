import numpy as np

# --- A dictionary containing galaxy-specific parameters taken from the literature ---
galaxies = {
    # EXAMPLE ENTRY
    "NGCexample": {"alpha_CO"           : 1.0,
                "distance"           : 100.0e6,
                "Ij"                 : np.array([10000, 20000, 30000]),
                "sigmaj"             : np.array([1.234, 2.345, 3.456]),
                "qj"                 : np.array([ 0.12,  0.23,  0.34]),
                "ML"                 : 2.0,
                "ra_center"          : 12.345678,
                "dec_center"         : -12.345678,
                "inclination"        : 50.0,
                "ellipse"            : {"x0_init": "center", "y0_init": "center", "sma_init": 0.25, "eps_init": 0.2, "pa_init": 45},
                "H2_mass"            : 10**9.45,
                "stellar_mass"       : 10**11.12,
                "log_stel_vel_disp"  : np.log10(250),
                "log_eff_radius"     : np.log10(10.0),
                "log_SFR"            : -0.50,
                "log_mu_star"        : 9.50,
                },

}