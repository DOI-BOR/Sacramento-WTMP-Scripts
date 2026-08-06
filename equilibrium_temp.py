"""
Equilibrium Temperature Utilities

This module implements shortwave/longwave radiation terms and iterative
root-finding routines (Newton–Raphson and bisection) to estimate the
water surface equilibrium temperature, adapted for Jython/Python use in
HEC-WAT / WTMP workflows.

"""

import datetime as dt   # Python datetime utilities for timestamps and arithmetic
import math             # Standard math library for trigonometric and exponential functions

import tz_offset        # Module providing timezone/day offset used in DOY calculations
reload(tz_offset)       # Reload to ensure the latest `tz_offset.days` value is used

# TODO: hardcoded lat/lon sould be variable

zero_K = 273.15                     # Celsius→Kelvin offset [K]
Tk = zero_K                         # Base Kelvin term used inside expansions (preserved as-is)
stefan_boltz_const = 5.67e-8        # Stefan–Boltzmann constant [W m^-2 K^-4]
lw_emissivity_water = 0.97          # Longwave emissivity of water [-]
reference_density = 1000.           # Reference water density [kg/m^3]
wind_a = 2.8e-9                     # Wind function coefficient (latent flux) [SI]
wind_b = 1.2e-9                     # Wind function coefficient (latent flux) [SI]
Magnus_a = 6.1078                   # Magnus formula coefficient (vapor pressure) [hPa]
Magnus_b = 17.27                    # Magnus formula coefficient [-]
Magnus_c = 237.7                    # Magnus formula coefficient [°C]
Cb = 0.61                           # Bulk transfer coefficient for sensible heat term [-]


def bisection_solve(fTe_min, a, b, H1, uw, Ta_in, td, sat_vp_ta, tol=1e-2, max_iter=1000):
    """
    Find a root using the bisection method for the function `fTe_min`.

    Parameters
    ----------
    fTe_min : callable
        Function of Te and additional parameters returning `Te - fTe(Te, ...)`.
    a : float
        Lower bound of the search interval.
    b : float
        Upper bound of the search interval.
    H1 : float
        Combined radiative heat flux term (shortwave + downwelling longwave) [W/m^2].
    uw : float
        Wind speed [m/s].
    Ta_in : float
        Air temperature [°C].
    td : float
        Dewpoint temperature [°C].
    sat_vp_ta : float
        Saturation vapor pressure at air temperature [hPa].
    tol : float, optional
        Convergence tolerance on interval half-width.
    max_iter : int, optional
        Maximum number of iterations before raising an error.

    Returns
    -------
    float
        Root estimate (equilibrium temperature) within `[a, b]` under tolerance.

    Notes
    -----
    - Guards against zero division when `(c - td) == 0` during iteration.

    Raises
    ------
    ValueError
        If maximum iterations are exceeded.
    """
    
    fa = fTe_min(a, H1, uw, Ta_in, td, sat_vp_ta)  # Evaluate at lower bound
    fb = fTe_min(b, H1, uw, Ta_in, td, sat_vp_ta)  # Evaluate at upper bound

    # Loop and solve until convergence
    iter_count = 0  # iteration counter
    while (b - a) / 2.0 > tol:
        # midpoint
        c = (a + b) / 2.0  
        
        # protect against divide-by-zero during iteration
        if (c - td) == 0.0:
            # nudge midpoint slightly
            c += tol  
        
        # evaluate at midpoint
        fc = fTe_min(c, H1, uw, Ta_in, td, sat_vp_ta)  
        
        # Check if tolerance is met
        if fc == 0 or (b - a) / 2.0 < tol:
            return c  # Root found within tolerance

        # Increment the iteration counter
        iter_count += 1
        if iter_count > max_iter:
            raise ValueError("bisection_solve: Maximum iterations exceeded")  # fail safe

        # Narrow interval based on sign change
        if fa * fc < 0.0:
            b = c
            fb = fc
        
        else:
            a = c
            fa = fc

    # Returning the middle of the final interval
    return (a + b) / 2.0  


def newton_raphson_solve(f, df, x0, H1, uw, Ta_in, Td, sat_vp_ta, tol=1e-2, max_iter=1000):
    """
    Newton–Raphson root finder for equilibrium temperature.

    Parameters
    ----------
    f : callable
        Function whose root is sought, signature `(Te, H1, uw, Ta_in, Td, sat_vp_ta)`.
    df : callable
        Derivative of `f` with respect to Te, same signature as `f`.
    x0 : float
        Initial guess for Te [°C].
    H1 : float
        Combined radiative heat flux term (shortwave + downwelling longwave) [W/m^2].
    uw : float
        Wind speed [m/s].
    Ta_in : float
        Air temperature [°C].
    Td : float
        Dewpoint temperature [°C].
    sat_vp_ta : float
        Saturation vapor pressure at air temperature [hPa].
    tol : float, optional
        Convergence tolerance on |x_new - x_old|.
    max_iter : int, optional
        Maximum iterations.

    Returns
    -------
    float
        Root estimate (equilibrium temperature) under tolerance.

    Raises
    ------
    ValueError
        If derivative magnitude is near zero or maximum iterations are exceeded.
    """
    
    # iteration counter
    iter_count = 0          
    
    # current iterate
    x_old = x0              
    
    # Iterate until the maximum number of iterations is met
    while iter_count < max_iter:
        # function value
        fx_old = f(x_old, H1, uw, Ta_in, Td, sat_vp_ta)   
        
        # derivative value
        dfx_old = df(x_old, H1, uw, Ta_in, Td, sat_vp_ta) 

        # Check if tolerance is exceeded
        if abs(dfx_old) < tol:
            raise ValueError("Derivative near zero")  # halt to avoid blow-up

        # Newton update
        x_new = x_old - fx_old / dfx_old  
        
        # Check for root value
        if abs(x_new - x_old) < tol:
            return x_new  # Root found

        # Advance the iteration and the counter
        x_old = x_new     
        iter_count += 1

    # Raise an error as this point should only be reached due to nonconvergence
    raise ValueError("newton_raphson_solve: Maximum iterations exceeded")


def get_decimal_day_of_year(tin):
    """
    Compute decimal day-of-year including timezone/day offset.

    Parameters
    ----------
    tin : datetime.datetime
        Timestamp (assumed local).

    Returns
    -------
    float
        Decimal DOY (1-based), adjusted by `tz_offset.days`.

    Notes
    -----
    The offset `tz_offset.days` is added directly to the DOY to reflect
    application-specific timezone shifts. Behavior is preserved as-is.
    """
    
    # First-of-year at midnight
    foy = dt.datetime(tin.year, 1, 1, 0, 0, 0)                      
    
    # Convert to days + offset
    ddoy = (tin - foy).total_seconds() / 86400. + 1. + tz_offset.days  
    
    # Return date to the calling function
    return ddoy


def solar_alt_angle(tin):
    """
    Calculate solar altitude angle (degrees) for a timestamp with fixed lat/lon.

    Parameters
    ----------
    tin : datetime.datetime
        Timestamp for which to compute solar altitude.

    Returns
    -------
    float
        Solar altitude angle in degrees.

    Notes
    -----
    Latitude, longitude, and time-zone longitude are hardcoded:
    """
    # todo: address the hard coded values here
    
    # Set the location information
    latitude = 38.1
    longitude = 121.8
    time_zone_longitude = 120.

    # Calculate the offset values
    delta = (longitude - time_zone_longitude) / 15.  # time-zone hour offset
    doy = get_decimal_day_of_year(tin)               # decimal day-of-year

    # Calculate the solar angle
    solar_decl = 0.4092 * math.cos(2. * math.pi / 365. * (172. - doy))  # declination [rad]
    
    t1 = math.sin(solar_decl) * math.sin(latitude * math.pi / 180.)      # component 1
    t2 = math.cos(solar_decl) * math.cos(latitude * math.pi / 180.)      # component 2
    
    hr = (doy - math.floor(doy)) * 24.                                     # local hour
    hr_angle = math.pi * (hr - delta - 12.) / 12.                          # hour angle [rad]
    
    sin_solar_alt = t1 + t2 * math.cos(hr_angle)                           # sine of altitude
    solar_alt_angle = math.asin(sin_solar_alt) * 180. / math.pi            # altitude [deg]
    
    # Return the solar angle to the calling function
    return solar_alt_angle


def sw_water_reflectance(tin, cloud_in):
    """
    Estimate shortwave water surface reflectance as a function of solar altitude and cloud cover.

    Parameters
    ----------
    tin : datetime.datetime
        Timestamp used to compute solar altitude.
    cloud_in : float
        Cloud fraction [0, 1].

    Returns
    -------
    float
        Water surface reflectance fraction [0, 1].
    """
    
    # Calculate the solare angle
    alt_angle = solar_alt_angle(tin)  # solar altitude
    
    # Adjust based on the angle magnitude
    if alt_angle <= 0.:
        return 1.  # at/below horizon: full reflection
    
    else:
        # Piecewise coefficients based on cloud cover regime
        if cloud_in > 0.95:
            acoef = 0.33
            bcoef = -0.45
    
        elif cloud_in > 0.55:
            acoef = 0.95
            bcoef = -0.75
        
        elif cloud_in > 0.05:
            acoef = 2.20
            bcoef = -0.97
        
        else:
            acoef = 1.18
            bcoef = -0.77

        # Empirical power-law reflectance limited to 1.0
        water_sfc_reflectance = acoef * alt_angle ** bcoef
        water_sfc_reflectance = min(water_sfc_reflectance, 1.)
        
        # Return to the calling function
        return water_sfc_reflectance


def heat_flux_surface_longwave_down(air_temp_in, cloud_in):
    """
    Compute downwelling longwave heat flux at the water surface.

    Parameters
    ----------
    air_temp_in : float
        Air temperature [°C].
    cloud_in : float
        Cloud fraction [0, 1].

    Returns
    -------
    float
        Downward longwave heat flux to surface [W/m^2].
    """
    
    # Set coefficient values
    sfc_reflect_fract = 0.03               # surface reflectance fraction for LW
    cloud_cover_fract_coef = 0.17          # cloud enhancement coefficient
    con1 = 0.937e-5                        # emissivity scaling coefficient
    
    # Convert to Kelvin
    Tkel = air_temp_in + zero_K            
    
    # Calculate the effective emissivity
    emissivity_air = con1 * (1. + cloud_cover_fract_coef * cloud_in ** 2) * Tkel ** 2 

    # Calculate teh downward long wave flux 
    hlwd = emissivity_air * stefan_boltz_const * (1. - sfc_reflect_fract) * Tkel ** 4 
    
    # Return the long way flux to the calling function
    return hlwd


def latent_heat_vaporization(temp_in):
    """
    Compute latent heat of vaporization as a function of temperature.

    Parameters
    ----------
    temp_in : float
        Temperature [°C].

    Returns
    -------
    float
        Latent heat of vaporization [J/kg].
    """
    
    return 1000. * (2499. - 2.36 * temp_in)  # linear approximation (°C dependence)


def sat_water_vapor_pres(temp_in):
    """
    Saturation water vapor pressure via Magnus–Tetens formula.

    Parameters
    ----------
    temp_in : float
        Temperature [°C].

    Returns
    -------
    float
        Saturation vapor pressure [hPa].
    """
    
    # Magnes-Tetens formula
    return Magnus_a * math.exp(Magnus_b * temp_in / (temp_in + Magnus_c))


def fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    """
    Equilibrium temperature residual function `fTe(Te)`.

    Parameters
    ----------
    Te : float
        Surface equilibrium temperature [°C].
    H1 : float
        Combined radiative term (shortwave + longwave down) [W/m^2].
    uw : float
        Wind speed [m/s].
    Ta_in : float
        Air temperature [°C].
    Td : float
        Dewpoint temperature [°C].
    sat_vp_ta : float
        Saturation vapor pressure at air temperature [hPa].

    Returns
    -------
    float
        Residual value for equilibrium condition at `Te`.

    Notes
    -----
    This function includes linearized longwave emission terms and combined
    sensible/latent flux terms using bulk transfer coefficients.
    """
    '''Calculation of fTe with only Te-dependent terms'''
    
    Lw = latent_heat_vaporization(Te)                     # latent heat (J/kg)
    sat_vp_te = sat_water_vapor_pres(Te)                  # saturation VP at Te
    beta = (sat_vp_te - sat_vp_ta) / (Te - Td)            # humidity slope term
    
    return H1 - lw_emissivity_water * stefan_boltz_const * (Tk ** 4 + 4. * Tk ** 3 * Te + 6. * Tk ** 2 * Te ** 2) + \
        reference_density * Lw * (wind_a + wind_b * uw) * ((Cb * Ta_in + beta * Td) - (Cb + beta) * Te)


def dfTe_dTe(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    """
    Derivative of `fTe` with respect to Te.

    Parameters
    ----------
    Te : float
        Surface equilibrium temperature [°C].
    H1 : float
        Combined radiative term (shortwave + longwave down) [W/m^2].
    uw : float
        Wind speed [m/s].
    Ta_in : float
        Air temperature [°C].
    Td : float
        Dewpoint temperature [°C].
    sat_vp_ta : float
        Saturation vapor pressure at air temperature [hPa].

    Returns
    -------
    float
        d(fTe)/dTe evaluated at `Te`.


    """


    Lw = 1000. * (2499. - 2.36 * Te)
    sat_vp_te = Magnus_a * math.exp(Magnus_b * Te / (Te + Magnus_c))

    # Derivatives
    dLw_dTe = -1000. * 2.36
    dsat_vp_te_dTe = (Magnus_a * Magnus_b * math.exp(Magnus_b * Te / (Te + Magnus_c))) / (Te + Magnus_c) - \
                     (Magnus_a * Magnus_b * Te * math.exp(Magnus_b * Te / (Te + Magnus_c))) / (Te + Magnus_c) ** 2
    dbeta_dTe = (dsat_vp_te_dTe * (Te - Td) - (sat_vp_te - sat_vp_ta)) / (Te - Td) ** 2

    dTe_term = - 4.0 * lw_emissivity_water * stefan_boltz_const * Tk ** 3 - 12.0 * lw_emissivity_water * stefan_boltz_const * Tk ** 2 * Te
    Lw_term = reference_density * (wind_a + wind_b * uw) * ((Cb * Ta_in + (sat_vp_te - sat_vp_ta) / (Te - Td) * Td) - (
            Cb + (sat_vp_te - sat_vp_ta) / (Te - Td)) * Te)
    beta_term = reference_density * Lw * (wind_a + wind_b * uw) * (-(Cb + (sat_vp_te - sat_vp_ta) / (Te - Td)))

    return dTe_term + Lw_term + beta_term + reference_density * dLw_dTe * (wind_a + wind_b * uw) * (
            (Cb * Ta_in + (sat_vp_te - sat_vp_ta) / (Te - Td) * Td) - (
            Cb + (sat_vp_te - sat_vp_ta) / (Te - Td)) * Te) + reference_density * Lw * (
            wind_a + wind_b * uw) * dbeta_dTe * Td


def fTe_min(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    """
    Helper returning `Te - fTe(Te, ...)` used by bisection search.

    Parameters
    ----------
    Te : float
        Trial equilibrium temperature [°C].
    H1, uw, Ta_in, Td, sat_vp_ta : see `fTe`

    Returns
    -------
    float
        `Te - fTe(Te, ...)` value.
    """
    
    return Te - fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta)


def fTe_abs(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    """
    Absolute residual magnitude used by numerical minimizers.

    Parameters
    ----------
    Te : float
        Trial equilibrium temperature [°C].
    H1, uw, Ta_in, Td, sat_vp_ta : see `fTe`

    Returns
    -------
    float
        `abs(fTe(Te, ...))`.
    """
    
    return abs(fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta))


def equilibrium_temp(dtt1, at, cl, sr, ws, td, te_guess, type='nr'):
    """
    Compute equilibrium temperature using selected solver.

    Parameters
    ----------
    dtt1 : datetime.datetime
        Timestamp for solar/reflectance terms.
    at : float
        Air temperature [°C].
    cl : float
        Cloud fraction [0, 1].
    sr : float
        Shortwave irradiance [W/m^2].
    ws : float
        Wind speed [m/s].
    td : float
        Dewpoint temperature [°C].
    te_guess : float
        Initial guess for Newton–Raphson (or SciPy minimize).
    type : {'nr', 'bs', 'scipy'}, optional
        Solver type: Newton–Raphson ('nr'), bisection ('bs'), or SciPy ('scipy').

    Returns
    -------
    float
        Equilibrium temperature [°C].

    Raises
    ------
    ValueError
        If `type` is not recognized.

    Notes
    -----
    - The SciPy path is preserved (commented rationale) but not usable in Jython.
    - Shortwave reflectance depends on solar altitude and cloud fraction.
    """
    
    reflectance = sw_water_reflectance(dtt1, cl)                    # SW reflectance
    Hsw = sr * (1. - reflectance)                                   # net shortwave to water
    HH = heat_flux_surface_longwave_down(at, cl)                    # downwelling longwave
    H1 = Hsw + HH                                                   # combined radiative term
    sat_vp_ta = sat_water_vapor_pres(at)                            # saturation VP at air T

    # Take action based on the specified type
    if type == 'nr':
        return newton_raphson_solve(fTe, dfTe_dTe, te_guess, H1, ws, at, td, sat_vp_ta, tol=1.0e-2, max_iter=1000)
    
    elif type == 'bs':
        return bisection_solve(fTe_min, -30.0, 55.0, H1, ws, at, td, sat_vp_ta, tol=1.0e-2, max_iter=1000)
    
    elif type == 'scipy':
       # Original equilibrium temp calc using scipy.optimize/minimize.  Can't use scipy in jython
       args = (H1, ws, at, td, sat_vp_ta)
       res = minimize(fTe_abs, x0=te_guess, args=args)
       return res.x
    
    else:
        raise ValueError('equilibrium_temp() - calculation type not known:', type)


def calc_equilibrium_temp(dtt, at, cl, sr, td, ws):
    """
    Vectorized wrapper to compute equilibrium temperature through time.

    Parameters
    ----------
    dtt : list of datetime.datetime
        Time series of timestamps.
    at : list of float
        Air temperature series [°C].
    cl : list of float
        Cloud fraction series [0, 1].
    sr : list of float
        Shortwave irradiance series [W/m^2].
    td : list of float
        Dewpoint temperature series [°C].
    ws : list of float
        Wind speed series [m/s].

    Returns
    -------
    list of float
        Equilibrium temperature series [°C].

    Notes
    -----
    Preserves behavior of trying Newton–Raphson, then using
    bisection as a fallback and sanity check (replace if |difference| > 1°C).
    """

    # Define a holder for the data
    nt = len(dtt)         # number of timesteps
    te = []               # output equilibrium temperature list

    # Loop and solve for each timestep
    for j in range(nt):

        # Handle the guess based on the timestep
        if j == 0:
            x0 = at[j]    # initial guess from air temperature
        else:
            x0 = te[j - 1]  # use last solution as next initial guess

        # The Newton-Raphson solver sometimes fails (when the derivitive goes big/non-linear) or produces a really bad
        # value.  We calulcate a secondary solution (which is sometimes inaccurate be over a degree, but usually is
        # within 0.2 degrees of the solution using scipy.optimize.minimize) and use it in those cases.
        te.append(-999)   # placeholder value
        te_bs = equilibrium_temp(dtt[j], at[j], cl[j], sr[j], ws[j], td[j], x0, type='bs')  # bisection estimate
        try:
            te[j] = equilibrium_temp(dtt[j], at[j], cl[j], sr[j], ws[j], td[j], x0, type='nr')  # Newton–Raphson
        except:
            # print(j, ' Newton-Raphson failure (convergence).  Using bisection-solution equilibrium temp')
            te[j] = te_bs  # fallback to bisection

        if abs(te[j] - te_bs) > 1.0:
            te[j] = te_bs  # replace with bisection when discrepancy is large

    # Return the equilibrium temperature to the calling function
    return te