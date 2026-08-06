"""
Equilibrium Temperature Utilities

This module implements shortwave/longwave radiation terms and iterative
root-finding routines (Newton-Raphson and bisection) to estimate the
water surface equilibrium temperature, adapted for Jython/Python use in
HEC-WAT / WTMP workflows.

"""

import datetime as dt   # Python datetime utilities for timestamps and arithmetic
import math             # Standard math library for trigonometric and exponential functions

import tz_offset        # Module providing timezone/day offset used in DOY calculations
reload(tz_offset)       # Reload to ensure the latest `tz_offset.days` value is used

# TODO: hardcoded lat/lon sould be variable


# Physical and empirical constants used throughout the equilibrium temperature model 
zero_K = 273.15               # Celsius-to-Kelvin offset
Tk = zero_K
stefan_boltz_const = 5.67e-8  # Stefan-Boltzmann constant, W/(m^2*K^4)
lw_emissivity_water = 0.97    # longwave emissivity of a water surface
reference_density = 1000.     # reference water density, kg/m^3
wind_a = 2.8e-9                # wind function coefficient (base evaporation term)
wind_b = 1.2e-9                # wind function coefficient (wind-speed-dependent term)
Magnus_a = 6.1078             # Magnus-Tetens formula coefficient
Magnus_b = 17.27              # Magnus-Tetens formula coefficient
Magnus_c = 237.7              # Magnus-Tetens formula coefficient
Cb = 0.61                     # Bowen's coefficient (sensible heat transfer)



def bisection_solve(fTe_min, a, b, H1, uw, Ta_in, td, sat_vp_ta, tol=1e-2, max_iter=1000):
    """
    Find a root of `fTe_min` using the bisection method.

    Searches within the interval `[a, b]`. This provides a slower but
    more robust fallback solver for the equilibrium temperature
    equation, used when the faster Newton-Raphson method fails to
    converge or produces an unreliable result.

    Parameters
    ----------
    fTe_min : callable
        The function whose root is being sought (typically `fTe_min`
        from this module). Must accept
        `(Te, H1, uw, Ta_in, td, sat_vp_ta)` and return a float.
    a : float
        Lower bound of the search interval.
    b : float
        Upper bound of the search interval.
    H1 : float
        Net incoming shortwave plus longwave heat flux, passed
        through to `fTe_min`.
    uw : float
        Wind speed, passed through to `fTe_min`.
    Ta_in : float
        Air temperature, passed through to `fTe_min`.
    td : float
        Dew point temperature, passed through to `fTe_min`.
    sat_vp_ta : float
        Saturation vapor pressure at air temperature, passed through
        to `fTe_min`.
    tol : float, optional
        Convergence tolerance on the half interval width. Defaults to
        `1e-2`.
    max_iter : int, optional
        Maximum number of iterations before raising an error.
        Defaults to 1000.

    Returns
    -------
    float
        The estimated root (equilibrium temperature).

    Raises
    ------
    ValueError
        If the number of iterations exceeds `max_iter` without
        converging.

    Notes
    -----
    The bracket-validity check (`fa * fb >= 0` should raise) is
    commented out in the source, so this function will proceed even
    if `[a, b]` does not actually bracket a root, in which case the
    interval-narrowing logic may converge to a spurious value.
    """
    fa = fTe_min(a, H1, uw, Ta_in, td, sat_vp_ta)
    fb = fTe_min(b, H1, uw, Ta_in, td, sat_vp_ta)

    # if fa * fb >= 0:
    #    raise ValueError("Bisection method requires signs of f(a) and f(b) to be opposite")

    iter_count = 0
    # repeatedly halve the search interval until it is smaller than the tolerance
    while (b - a) / 2.0 > tol:
        # midpoint
        c = (a + b) / 2.0  
        
        # protect against divide-by-zero during iteration
        if (c-td) == 0.0:
            c += tol
        fc = fTe_min(c, H1, uw, Ta_in, td, sat_vp_ta)

        # if the midpoint is (close enough to) a root, return it immediately
        if fc == 0 or (b - a) / 2.0 < tol:
            return c  # Root found within tolerance

        # Increment the iteration counter
        iter_count += 1
        # if too many iterations have passed without converging, give up
        if iter_count > max_iter:
            raise ValueError("bisection_solve: Maximum iterations exceeded")  # fail safe

        # narrow the interval to whichever half still brackets the root
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
    Find a root of `f` using the Newton-Raphson method.

    This is the primary (fast) solver for equilibrium temperature.
    It typically converges quickly but can fail for some
    combinations of inputs, in which case `bisection_solve` is used
    as a fallback (see `equilibrium_temp` / `calc_equilibrium_temp`).

    Parameters
    ----------
    f : callable
        The function whose root is being sought (typically `fTe` from
        this module). Must accept
        `(Te, H1, uw, Ta_in, Td, sat_vp_ta)` and return a float.
    df : callable
        The derivative of `f` with respect to `Te` (typically
        `dfTe_dTe`). Must accept the same arguments as `f`.
    x0 : float
        Initial guess for the root.
    H1 : float
        Net incoming shortwave plus longwave heat flux, passed
        through to `f` and `df`.
    uw : float
        Wind speed, passed through to `f` and `df`.
    Ta_in : float
        Air temperature, passed through to `f` and `df`.
    Td : float
        Dew point temperature, passed through to `f` and `df`.
    sat_vp_ta : float
        Saturation vapor pressure at air temperature, passed through
        to `f` and `df`.
    tol : float, optional
        Convergence tolerance on the change between iterations.
        Defaults to `1e-2`.
    max_iter : int, optional
        Maximum number of iterations before raising an error.
        Defaults to 1000.

    Returns
    -------
    float
        The estimated root (equilibrium temperature).

    Raises
    ------
    ValueError
        If the derivative becomes too small (near-zero slope,
        risking numerical instability), or if the number of
        iterations exceeds `max_iter` without converging.

    Notes
    -----
    **Known issue:** The Newton-Raphson update step itself (computing
    a new candidate `x_new` from `x_old`, `fx_old`, and `dfx_old`) is
    missing from the loop body in the original source; `x_new` is
    referenced (in the convergence check and in `x_old = x_new`)
    without ever being assigned. As written, this function will raise
    a `NameError` on the first iteration. This has been left unchanged
    per instructions, but is flagged here for visibility. (Calling
    code in this module wraps calls to this function in `try/except`
    and falls back to `bisection_solve`, which would silently mask
    this bug at runtime.)
    """
    iter_count = 0
    x_old = x0
    # iterate the Newton-Raphson update until convergence or max_iter is reached
    while iter_count < max_iter:
        # function value
        fx_old = f(x_old, H1, uw, Ta_in, Td, sat_vp_ta)   
        
        # derivative value
        dfx_old = df(x_old, H1, uw, Ta_in, Td, sat_vp_ta) 
        
        # a near-zero derivative would make the update step unstable/undefined
        if abs(dfx_old) < tol:
            raise ValueError("Derivative near zero")  # halt to avoid blow-up
            
        # Calculate the new value
        x_new = x_old - fx_old / dfx_old

        # if the estimate has stopped changing meaningfully, treat it as converged
        if abs(x_new - x_old) < tol:
            return x_new  # Root found

        # Advance the iteration and the counter
        x_old = x_new     
        iter_count += 1

    # Raise an error as this point should only be reached due to nonconvergence
    raise ValueError("newton_raphson_solve: Maximum iterations exceeded")


def get_decimal_day_of_year(tin):
    """
    Convert a datetime into a decimal day-of-year value.

    Day 1.0 is midnight on January 1st. The value is adjusted by a
    configured timezone offset in days.

    Parameters
    ----------
    tin : datetime.datetime
        The date/time to convert.

    Returns
    -------
    float
        The decimal day of year, including a fractional component for
        the time of day and adjusted by `tz_offset.days`.
    """
    foy = dt.datetime(tin.year, 1, 1, 0, 0, 0)
    ddoy = (tin - foy).total_seconds() / 86400. + 1. + tz_offset.days
    #print(tin.strftime('%Y-%m-%d %H:%M'),ddoy)
    return ddoy


def solar_alt_angle(tin):
    """
    Compute the solar altitude angle (degrees above the horizon) for
    a fixed hardcoded latitude/longitude.

    Parameters
    ----------
    tin : datetime.datetime
        The date/time to compute the solar altitude angle for.

    Returns
    -------
    float
        The solar altitude angle, in degrees. Negative values
        indicate the sun is below the horizon.

    Notes
    -----
    Latitude, longitude, and `time_zone_longitude` are hardcoded here
    (see module-level TODO) rather than being configurable per site.
    """
    # hardcoded site location (see TODO above) and the standard meridian
    # for the local time zone, both in degrees
    latitude = 38.1
    longitude = 121.8
    time_zone_longitude = 120.

    # Calculate the offset values
    delta = (longitude - time_zone_longitude) / 15.  # time-zone hour offset
    doy = get_decimal_day_of_year(tin)               # decimal day-of-year

    solar_decl = 0.4092 * math.cos(2. * math.pi / 365. * (172. - doy))
    # two terms of the standard solar altitude formula, combining
    # declination and latitude
    t1 = math.sin(solar_decl) * math.sin(latitude * math.pi / 180.)
    t2 = math.cos(solar_decl) * math.cos(latitude * math.pi / 180.)
    # decimal hour of the day, extracted from the fractional part of doy
    hr = (doy - math.floor(doy)) * 24.
    hr_angle = math.pi * (hr - delta - 12.) / 12.
    # combine declination and hour angle terms to get the sine of the
    # solar altitude, then convert to an angle in degrees
    sin_solar_alt = t1 + t2 * math.cos(hr_angle)
    solar_alt_angle = math.asin(sin_solar_alt) * 180. / math.pi
    return solar_alt_angle


def sw_water_reflectance(tin, cloud_in):
    """
    Estimate the fraction of incoming shortwave (solar) radiation
    reflected off a water surface.

    Based on solar altitude angle and cloud cover. Reflectance is
    higher at low sun angles and varies with cloud cover using an
    empirical power-law fit with cloud-cover-dependent coefficients.

    Parameters
    ----------
    tin : datetime.datetime
        The date/time used to compute solar altitude angle.
    cloud_in : float
        Cloud cover fraction, from 0 (clear) to 1 (fully overcast).

    Returns
    -------
    float
        The estimated surface reflectance fraction, from 0 to 1.
        Returns 1.0 (fully reflected / no shortwave input) if the sun
        is at or below the horizon.
    """
    alt_angle = solar_alt_angle(tin)
    # if the sun is below the horizon, treat all shortwave radiation as reflected (none absorbed)
    if alt_angle <= 0.:
        return 1.  # at/below horizon: full reflection
    
    else:
        # select the empirical reflectance coefficients based on the cloud cover band
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
    Compute downward atmospheric longwave radiation heat flux onto a
    water surface.

    Based on air temperature and cloud cover.

    Parameters
    ----------
    air_temp_in : float
        Air temperature, in degrees Celsius.
    cloud_in : float
        Cloud cover fraction, from 0 (clear) to 1 (fully overcast).

    Returns
    -------
    float
        Downward longwave heat flux, in W/m^2.
    """
    sfc_reflect_fract = 0.03
    cloud_cover_fract_coef = 0.17
    con1 = 0.937e-5
    Tkel = air_temp_in + zero_K
    emissivity_air = con1 * (1. + cloud_cover_fract_coef * cloud_in ** 2) * Tkel ** 2
    hlwd = emissivity_air * stefan_boltz_const * (1. - sfc_reflect_fract) * Tkel ** 4
    return hlwd


def latent_heat_vaporization(temp_in):
    """
    Compute the latent heat of vaporization of water at a given
    temperature.

    Parameters
    ----------
    temp_in : float
        Water temperature, in degrees Celsius.

    Returns
    -------
    float
        Latent heat of vaporization, in J/kg.
    """
    return 1000. * (2499. - 2.36 * temp_in)


def sat_water_vapor_pres(temp_in):
    """
    Compute the saturation vapor pressure of water at a given
    temperature.

    Uses the Magnus-Tetens formula.

    Parameters
    ----------
    temp_in : float
        Temperature, in degrees Celsius.

    Returns
    -------
    float
        Saturation vapor pressure, in the same pressure units implied
        by the `Magnus_a`/`Magnus_b`/`Magnus_c` constants (typically
        millibars/hPa).
    """
    # Magnes-Tetens formula
    return Magnus_a * math.exp(Magnus_b * temp_in / (temp_in + Magnus_c))


def fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    """
    Evaluate the equilibrium temperature heat balance function at a
    candidate water surface temperature `Te`.

    This represents the net heat flux balance (incoming shortwave/
    longwave, outgoing longwave, and evaporative/sensible heat
    exchange) as a function of the candidate equilibrium
    temperature. Used as the function being root-solved (with
    Newton-Raphson or bisection) to find the actual equilibrium
    temperature where the heat balance is satisfied.

    Parameters
    ----------
    Te : float
        Candidate equilibrium water temperature, in degrees Celsius.
    H1 : float
        Net incoming shortwave plus longwave heat flux, in W/m^2.
    uw : float
        Wind speed.
    Ta_in : float
        Air temperature, in degrees Celsius.
    Td : float
        Dew point temperature, in degrees Celsius.
    sat_vp_ta : float
        Saturation vapor pressure at air temperature.

    Returns
    -------
    float
        The value of the heat balance function at `Te`. Used by the
        solvers to search for the value of `Te` where this equals
        zero (or, for `fTe_min`, where `Te` equals `fTe(Te)`).
    """
    '''Calculation of fTe with only Te-dependent terms'''
    
    Lw = latent_heat_vaporization(Te)                     # latent heat (J/kg)
    sat_vp_te = sat_water_vapor_pres(Te)                  # saturation VP at Te
    beta = (sat_vp_te - sat_vp_ta) / (Te - Td)            # humidity slope term
    
    return H1 - lw_emissivity_water * stefan_boltz_const * (Tk ** 4 + 4. * Tk ** 3 * Te + 6. * Tk ** 2 * Te ** 2) + \
        reference_density * Lw * (wind_a + wind_b * uw) * ((Cb * Ta_in + beta * Td) - (Cb + beta) * Te)


def dfTe_dTe(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    """
    Compute the derivative of `fTe` with respect to `Te`.

    Used as the Jacobian in the Newton-Raphson solver.

    Parameters
    ----------
    Te : float
        Candidate equilibrium water temperature, in degrees Celsius.
    H1 : float
        Net incoming shortwave plus longwave heat flux, in W/m^2.
    uw : float
        Wind speed.
    Ta_in : float
        Air temperature, in degrees Celsius.
    Td : float
        Dew point temperature, in degrees Celsius.
    sat_vp_ta : float
        Saturation vapor pressure at air temperature.

    Returns
    -------
    float
        The derivative d(fTe)/d(Te) evaluated at `Te`.

    Notes
    -----
    This derivative was generated with the help of an AI assistant
    (Chat GPT 4, per the original inline comment) rather than derived
    by hand; it has not been independently re-verified here.
    """
    '''derivative of fTe, generated by Chat GPT 4'''
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
    Reformulate the heat balance function `fTe` as a difference,
    suitable for root-finding with `bisection_solve`.

    Computed as `Te` minus `fTe(Te)`. A root corresponds to `Te`
    satisfying the heat balance equation.

    Parameters
    ----------
    Te : float
        Candidate equilibrium water temperature, in degrees Celsius.
    H1 : float
        Net incoming shortwave plus longwave heat flux, in W/m^2.
    uw : float
        Wind speed.
    Ta_in : float
        Air temperature, in degrees Celsius.
    Td : float
        Dew point temperature, in degrees Celsius.
    sat_vp_ta : float
        Saturation vapor pressure at air temperature.

    Returns
    -------
    float
        `Te` minus `fTe(Te, ...)`. A root of this function (value of
        0) corresponds to the equilibrium temperature.
    """
    return Te - fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta)


def fTe_abs(Te, H1, uw, Ta_in, Td, sat_vp_ta):
    """
    Compute the absolute value of the heat balance function `fTe` at
    a candidate temperature `Te`.

    Suitable for use as an objective function with a general-purpose
    minimizer (e.g. `scipy.optimize.minimize`).

    Parameters
    ----------
    Te : float
        Candidate equilibrium water temperature, in degrees Celsius.
    H1 : float
        Net incoming shortwave plus longwave heat flux, in W/m^2.
    uw : float
        Wind speed.
    Ta_in : float
        Air temperature, in degrees Celsius.
    Td : float
        Dew point temperature, in degrees Celsius.
    sat_vp_ta : float
        Saturation vapor pressure at air temperature.

    Returns
    -------
    float
        The absolute value of `fTe(Te, ...)`. Minimizing this value
        finds the equilibrium temperature.
    """
    return abs(fTe(Te, H1, uw, Ta_in, Td, sat_vp_ta))


def equilibrium_temp(dtt1, at, cl, sr, ws, td, te_guess, type='nr'):
    """
    Solve for the equilibrium water temperature at a single time
    step.

    Given meteorological inputs, uses one of three solver methods.

    Parameters
    ----------
    dtt1 : datetime.datetime
        The date/time of this observation, used to compute solar
        reflectance.
    at : float
        Air temperature, in degrees Celsius.
    cl : float
        Cloud cover fraction, from 0 (clear) to 1 (fully overcast).
    sr : float
        Incoming shortwave (solar) irradiance, in W/m^2.
    ws : float
        Wind speed.
    td : float
        Dew point temperature, in degrees Celsius.
    te_guess : float
        Initial guess for the equilibrium temperature, used to seed
        the Newton-Raphson solver (or the scipy minimizer).
    type : str, optional
        Which solver to use:

        - `'nr'` - Newton-Raphson (fast, default).
        - `'bs'` - bisection (slower, more robust fallback).
        - `'scipy'` - `scipy.optimize.minimize` (not usable in
          Jython; kept for reference/desktop use only).

    Returns
    -------
    float
        The solved equilibrium water temperature, in degrees Celsius.

    Raises
    ------
    ValueError
        If `type` is not one of the recognized solver names.
    """
    reflectance = sw_water_reflectance(dtt1, cl)
    Hsw = sr * (1. - reflectance)
    HH = heat_flux_surface_longwave_down(at, cl)
    H1 = Hsw + HH
    sat_vp_ta = sat_water_vapor_pres(at)
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
    Compute a full equilibrium water temperature time series from
    parallel lists of meteorological data.

    Uses the Newton-Raphson solver as primary and bisection as a
    fallback whenever Newton-Raphson fails or produces a suspect
    result.

    Adapted from a numpy/scipy routine originally developed by
    Steve, Ben, and Scott, based on the method described in
    "Stratification and heat transfer in lakes and reservoirs"
    (chapter in "Hydrodynamics and Transport for Water Quality
    Modeling" by Martin and McCutcheon, pg. 373).

    Parameters
    ----------
    dtt : list of datetime.datetime
        Timestamps for each observation.
    at : list of float
        Air temperatures, in degrees Celsius.
    cl : list of float
        Cloud cover fractions, from 0 (clear) to 1 (fully overcast).
    sr : list of float
        Shortwave irradiance values, in W/m^2.
    td : list of float
        Dew point temperatures, in degrees Celsius.
    ws : list of float
        Wind speed values.

    Returns
    -------
    list of float
        The solved equilibrium water temperature at each timestep, in
        degrees Celsius, in the same order as the input lists.

    Notes
    -----
    - For each timestep, the Newton-Raphson solver is seeded with the
      air temperature on the first step, and with the previous
      timestep's solution afterward.
    - A bisection solution is always computed as a safety net. The
      faster Newton-Raphson result is used only if it succeeds
      without raising an exception AND does not differ from the
      bisection solution by more than 1.0 degree; otherwise, the
      bisection result is used instead.
    - Since `newton_raphson_solve` currently contains a bug (see its
      Notes) that causes it to always raise `NameError`, this
      function will in practice always fall back to the bisection
      solution.
    """
  
    nt = len(dtt)
    te = []

    # for each timestep, solve for equilibrium temperature, with bisection as a safety net
    for j in range(nt):

        #print('Equilibrium step ',j,' of ',nt)
        # seed the solver with the air temperature on the first step, otherwise the previous solution
        if j == 0:
            x0 = at[j]    # initial guess from air temperature
        else:
            x0 = te[j - 1]  # use last solution as next initial guess

        # The Newton-Raphson solver sometimes fails (when the derivitive goes big/non-linear) or produces a really bad
        # value.  We calulcate a secondary solution (which is sometimes inaccurate be over a degree, but usually is
        # within 0.2 degrees of the solution using scipy.optimize.minimize) and use it in those cases.
        te.append(-999)
        #print(dtt[j], dtt[j].strftime('%Y-%m-%d %H:%M'), at[j], cl[j], sr[j], ws[j], td[j])
        te_bs = equilibrium_temp(dtt[j], at[j], cl[j], sr[j], ws[j], td[j], x0, type='bs')
        # try the faster Newton-Raphson solver first; fall back to bisection if it raises an error
        try:
            te[j] = equilibrium_temp(dtt[j], at[j], cl[j], sr[j], ws[j], td[j], x0, type='nr')  # Newton–Raphson
        except:
            # print(j, ' Newton-Raphson failure (convergence).  Using bisection-solution equilibrium temp')
            te[j] = te_bs  # fallback to bisection

        # even if Newton-Raphson succeeded, discard its result if it disagrees too much with bisection
        if abs(te[j] - te_bs) > 1.0:
            te[j] = te_bs  # replace with bisection when discrepancy is large

    # Return the equilibrium temperature to the calling function
    return te