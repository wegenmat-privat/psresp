import numpy as np

__all__ = [
    'psresp',
]


def _spectrum(x, slope):
    y = x**slope

    return y


def _timmerlc(slope, nt='None', dt='None', mean='None', sigma='None', seed='None'):
    if dt == 'None':
        dt = 1
    if nt == 'None':
        nt = 65536
    if mean == 'None':
        mean = 0
    if sigma == 'None':
        sigma = 1
    if seed == 'None':
        seed = 42

    simfreq = np.linspace(1, nt/2-1, num=nt/2, dtype='float64') / (dt*nt)
    simpsd = _spectrum(simfreq, slope)
    fac = np.sqrt(simpsd)

    pos_real = np.random.RandomState(seed).normal(size=int(nt/2))*fac
    pos_imag = np.random.RandomState(seed).normal(size=int(nt/2))*fac

    pos_imag[int(nt/2)-1] = 0

    if float(nt/2.) > int(nt/2):
        neg_real = pos_real[0:int(nt/2)][::-1]
        neg_imag = -pos_real[0:int(nt/2)][::-1]
    else:
        neg_real = pos_real[0:int(nt/2)-1][::-1]
        neg_imag = -pos_real[0:int(nt/2)-1][::-1]

    real = np.hstack((0., pos_real, neg_real))
    imag = np.hstack((0., pos_imag, neg_imag))

    arg = real + 1j * imag
    rate = np.fft.ifft(arg).real
    time = dt*np.linspace(0, nt-1, nt, dtype='float')

    avg = np.mean(rate)
    std = np.sqrt(np.var(rate))

    rate = (rate - avg) * sigma / std + mean

    return time, rate


def _periodogram(series):
    series = series - np.mean(series)
    revff = np.fft.ifft(series)
    periodogram_bar = np.zeros(len(revff))
    freqs = np.zeros(len(revff))
    for freq, ff in enumerate(revff):
        if freq == 0:
            continue
        freqs[freq] = freq
        periodogram_bar[freq] = (abs(ff)**2)

    return freqs, periodogram_bar


def _rebinlc(time, rate, dt):
    # check for finite values
    rate = rate[np.isfinite(time)]
    time = time[np.isfinite(time)]

    # rebin lc/psd to evenly spaced time binning, from rebinlc.pro
    t = time

    ts = t[0]
    num = int((t[-1] - ts) / dt + 1)

    minnum = 1

    tn = np.empty([num])
    rn = np.empty([num])
    numperbin = np.empty([num])

    tlimit = ts+dt*np.linspace(0, num, num+1)

    k = 0
    for i in range(num):
        tn[k] = tlimit[i]
        index = np.where(np.greater_equal(t, tlimit[i]) & (t < tlimit[i+1]))
        number = len(t[index])

        numperbin[k] = number
        if np.greater_equal(number, minnum):
            rn[k] = np.sum(rate[index]) / number
            k = k+1
        else:
            rn[k] = 0
            numperbin[k] = 0
            k = k+1

    if k == 0:
        print('Rebinned lightcurve would not contain any data')

    if k != num:
        tn = tn[0:k]
        rn = rn[0:k]
        numperbin = numperbin[0:k]

    tn = tn[numperbin != 0]
    rn = rn[numperbin != 0]

    return tn, rn


def _do_periodogram(y):
    y = y - np.mean(y)
    freq, psd = _periodogram(y)
    
    return freq, psd


def _chi2_obs(norm, obs_pds, avg_pds, pds_err):
    obs_pds = norm * obs_pds
    a = (avg_pds - obs_pds)**2.
    b = pds_err**2.
    chi_obs = np.sum(a / b)

    return chi_obs


def _compare(obs_pds, avg_pds, pds_err, allpds, number_simulations):
    from scipy import optimize
    norm0 = 1.
    chi_obs = optimize.minimize(_chi2_obs, norm0, args=(obs_pds, avg_pds, pds_err), method='SLSQP').fun
    
    chi_dist = np.empty([number_simulations])
    for n in range(number_simulations):
        a = (allpds[n, :] - avg_pds)**2
        b = pds_err**2
        chi = a / b
        chi_dist[n] = np.sum(chi)

    suf = 0
    for i in range(len(chi_dist)):
        if np.greater_equal(chi_dist[i], chi_obs):
            suf = suf + 1
    
    suf = suf / len(chi_dist)

    return suf


def _psresp_pro(t, y, dy, slopes, number_simulations, binning, oversampling, df):
    bin = binning
    date, rat, raterr = t, y, dy
    date = date - date[0]
    duration = np.max(date) - np.min(date)
    npoints = int(duration / bin) * number_simulations * oversampling
    params = -slopes
    lc_variance = np.var(rat) - np.var(raterr)
    # observed PDS calculation
    obs_nu, obs_pds = _do_periodogram(rat)
    obs_nu = np.log10(obs_nu)

    # normalisation
    obs_pds = (2.*duration) / (np.mean(rat)*np.mean(rat)*len(rat)*len(rat)) * obs_pds

    # rebin
    obs_freqs, obs_power = _rebinlc(obs_nu, obs_pds, dt=df)
    obs_power = np.log10(obs_power)    

    # create fake light curve
    faketime, fakerate = _timmerlc(params, nt=npoints, dt=bin/oversampling)

    # calculate level of Poisson noise
    factor = ((len(raterr) / (2.*duration)) - (1. / duration))
    p_noise = np.sum(raterr**2.) / (len(raterr) * factor)

    # calculate high frequency aliased power
    uplim = 1. / (2.*bin)
    lowlim = 1. / (2.*(bin/10))
    intfreq = np.empty([int((lowlim-uplim)/uplim)+2])
    for i in range(len(intfreq)):
        intfreq[i] = uplim*(i+1)
    intpds = _spectrum(intfreq, params)
    integral = np.trapz(intpds, x=intfreq)
    p_alias = integral / factor

    # long light curve is divided and resultant PDS are calculated
    allpds = np.empty([number_simulations, len(obs_freqs)])
    for j in range(number_simulations):
        # print('computing PDS ' + str(j+1))
        
        # indices for each segment
        lobin = int(duration * j / (bin/oversampling))
        hibin = int(duration * j / (bin/oversampling)) + int(duration/(bin/oversampling))

        # taken from appropriate section of light curve
        temptime = faketime[lobin:hibin]
        temprate = fakerate[lobin:hibin]

        # shift start time to zero
        temptime = temptime - temptime[0]

        # set bintime equal to original light curve time
        bintime = date
        binrate = np.interp(date, temptime, temprate)        
        
        # rescale simulated LC to the mean and variance of the original
        tempvar = np.sqrt(np.var(binrate))
        binrate = (binrate - np.mean(binrate)) * ((np.sqrt(lc_variance)) / tempvar) + np.mean(rat)

        # calculate PDS of simulated light curve
        sim_nu, sim_pds = _do_periodogram(binrate)
        sim_nu = np.log10(sim_nu)
        # sim_pds = sim_pds + p_noise + p_alias

        # normalisation
        sim_pds = (2.*(np.max(bintime)-np.min(bintime))) / (np.mean(binrate)
                                                            * np.mean(binrate)*len(binrate)*len(binrate)) * sim_pds

        # rebin simulated PDS in same manner as observed
        # logfreqs, logpower = binlogPSD(sim_nu, sim_pds, df)
        logfreqs, power = _rebinlc(sim_nu, sim_pds, dt=df)
        logpower = np.log10(power)

        # large array for later calculations of mean and rms error
        for k in range(len(logpower)):
            allpds[j, k] = logpower[k]
    
    avg_pds = np.empty([len(obs_freqs)])    
    pds_err = np.empty([len(avg_pds)])
    for i in range(len(avg_pds)):
        avg_pds[i] = np.mean(allpds[:, i])
    for i in range(len(pds_err)):
        pds_err[i] = np.sqrt(np.var(allpds[:, i]))

    return faketime, fakerate, obs_freqs, obs_power, bintime, binrate, avg_pds, pds_err, allpds


def psresp(t, y, dy, slopes, dt, df, percentile, oversampling=10, number_simulations=100):
    """
    Compute power spectral density of a light curve assuming an unbroken power law with the PSRESP method.

    The artificial light curves are generated using the algorithm by Timmer and Koenig (1995).
    For an introduction to the PSRESP method, see Uttley (2002).

    The function returns a results dictionary with the following content:

    - ``slope`` (`float`) -- Mean slope of the power law
    - ``slope_error`` (`float`) -- Error of the slope of the power law
    - ``suf`` (`~numpy.ndarray`) -- Success fraction for each model parameter
    - ``best_parameters`` (`~numpy.ndarray`) -- Parameters satisfying the significance criterion
    - ``statistics`` (`~numpy.ndarray`) -- Data used to calculate the mean slope and its error
      over a grid of ``dt`` and ``df``

        - slope with the highest success fraction
        - highest success fraction
        - slope of the lower full width at half maximum for the success fraction distribution
        - slope of the higher full width at half maximum for the success fraction distribution

    Parameters
    ----------
    t : `~numpy.ndarray`
        Time array of the light curve
    y : `~numpy.ndarray`
        Flux array of the light curve
    dy : `~numpy.ndarray`
        Flux error array of the light curve
    slopes : `~numpy.ndarray`
        slopes of the power law model
    dt : `~numpy.ndarray`
        bin length for the light curve in units of ``t``
    df : `~numpy.ndarray`
        bin factor for the logarithmic periodogram
    percentile : `float`
        percentile of the distribution of success fraction, `0 < significance < 1`
    oversampling : `int`
        oversampling factor of the simulated light curve, default is 10
    number_simulations : `int`
        number of simulations for each model parameter, default is 10

    Returns
    -------
    results : `~collections.OrderedDict`
        Results dictionary (see description above).

    References
    ----------
    .. [1] Timmer and Koenig (1995), "On generating power law noise",
       `Link <http://adsabs.harvard.edu/abs/1995A%26A...300..707T>`_
    .. [2] Uttley et al, "Measuring the broad-band power spectra of active galactic nuclei with RXTE",
       `Link <https://academic.oup.com/mnras/article/332/1/231/974626/Measuring-the-broad-band-power-spectra-of-active>`_
    """

    from scipy import interpolate
    t_ini, y_ini, dy_ini = t, y, dy
    suf = np.empty([len(slopes), len(dt), len(df)])
    statistics = np.empty([4, len(dt), len(df)])
    for b in range(len(dt)):
        # print('binning: ' + str(dt[b]))
        t, y = _rebinlc(t_ini, y_ini, dt[b])
        t, dy = _rebinlc(t_ini, dy_ini, dt[b])
        for f in range(len(df)):
            # print('df: ' + str(df[f]))
            for s in range(len(slopes)):
                # print('slope: ' + str(slopes[s]))

                # psresp
                faketime, fakerate, obs_freqs, obs_power, bintime, binrate, avg_pds, pds_err, allpds = \
                    _psresp_pro(t, y, dy, slopes[s], number_simulations, dt[b], oversampling, df[f])

                # do chi2
                suf[s, b, f] = _compare(obs_power, avg_pds, pds_err, allpds, number_simulations)

            # find best slope and estimate error
            best_slope = slopes[np.argmax(suf[:, b, f])]
            best_slope_suf = np.max(suf[:, b, f])

            slopes_fwhm = interpolate.UnivariateSpline(slopes, suf[:, b, f] - 0.5 * best_slope_suf, s=0).roots()
            if len(slopes_fwhm) == 0:
                slopes_fwhm = [np.nan]
            low_slopes = slopes_fwhm[0]
            high_slopes = slopes_fwhm[-1]
            if low_slopes == high_slopes:
                low_slopes = high_slopes = np.nan

            statistics[0, b, f] = best_slope
            statistics[1, b, f] = best_slope_suf
            statistics[2, b, f] = low_slopes
            statistics[3, b, f] = high_slopes

    test_significance = np.percentile(statistics[1, :, :], 100*percentile)
    statistics_test = (np.greater_equal(statistics[1, :, :], test_significance)) & \
                      (np.isfinite(statistics[2, :, :])) & \
                      (np.isfinite(statistics[3, :, :]))
    best_parameters = np.where(statistics_test == True)
    mean_slope = np.sum(statistics[0, :, :][statistics_test] * statistics[1, :, :][statistics_test]) \
        / (np.sum(statistics_test))
    mean_error = np.abs(np.min(statistics[2, :, :][statistics_test]) - np.max(statistics[3, :, :][statistics_test]))

    return dict(slope=mean_slope,
                slope_error=mean_error,
                suf=suf,
                best_parameters=best_parameters,
                statistics=statistics
                )
