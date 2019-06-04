#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

from numpy.testing import (run_module_suite, assert_allclose, assert_warns,
                           assert_almost_equal, assert_raises)
import numpy as np
import pywt


def ref_gaus(LB, UB, N, num):
    X = np.linspace(LB, UB, N)
    F0 = (2./np.pi)**(1./4.)*np.exp(-(X**2))
    if (num == 1):
        psi = -2.*X*F0
    elif (num == 2):
        psi = -2/(3**(1/2))*(-1 + 2*X**2)*F0
    elif (num == 3):
        psi = -4/(15**(1/2))*X*(3 - 2*X**2)*F0
    elif (num == 4):
        psi = 4/(105**(1/2))*(3 - 12*X**2 + 4*X**4)*F0
    elif (num == 5):
        psi = 8/(3*(105**(1/2)))*X*(-15 + 20*X**2 - 4*X**4)*F0
    elif (num == 6):
        psi = -8/(3*(1155**(1/2)))*(-15 + 90*X**2 - 60*X**4 + 8*X**6)*F0
    elif (num == 7):
        psi = -16/(3*(15015**(1/2)))*X*(105 - 210*X**2 + 84*X**4 - 8*X**6)*F0
    elif (num == 8):
        psi = 16/(45*(1001**(1/2)))*(105 - 840*X**2 + 840*X**4 -
                                     224*X**6 + 16*X**8)*F0
    return (psi, X)


def ref_cgau(LB, UB, N, num):
    X = np.linspace(LB, UB, N)
    F0 = np.exp(-X**2)
    F1 = np.exp(-1j*X)
    F2 = (F1*F0)/(np.exp(-1/2)*2**(1/2)*np.pi**(1/2))**(1/2)
    if (num == 1):
        psi = F2*(-1j - 2*X)*2**(1/2)
    elif (num == 2):
        psi = 1/3*F2*(-3 + 4j*X + 4*X**2)*6**(1/2)
    elif (num == 3):
        psi = 1/15*F2*(7j + 18*X - 12j*X**2 - 8*X**3)*30**(1/2)
    elif (num == 4):
        psi = 1/105*F2*(25 - 56j*X - 72*X**2 + 32j*X**3 + 16*X**4)*210**(1/2)
    elif (num == 5):
        psi = 1/315*F2*(-81j - 250*X + 280j*X**2 + 240*X**3 -
                        80j*X**4 - 32*X**5)*210**(1/2)
    elif (num == 6):
        psi = 1/3465*F2*(-331 + 972j*X + 1500*X**2 - 1120j*X**3 - 720*X**4 +
                         192j*X**5 + 64*X**6)*2310**(1/2)
    elif (num == 7):
        psi = 1/45045*F2*(
            1303j + 4634*X - 6804j*X**2 - 7000*X**3 + 3920j*X**4 + 2016*X**5 -
            448j*X**6 - 128*X**7)*30030**(1/2)
    elif (num == 8):
        psi = 1/45045*F2*(
            5937 - 20848j*X - 37072*X**2 + 36288j*X**3 + 28000*X**4 -
            12544j*X**5 - 5376*X**6 + 1024j*X**7 + 256*X**8)*2002**(1/2)

    psi = psi/np.real(np.sqrt(np.real(np.sum(psi*np.conj(psi)))*(X[1] - X[0])))
    return (psi, X)


def sinc2(x):
    y = np.ones_like(x)
    k = np.where(x)[0]
    y[k] = np.sin(np.pi*x[k])/(np.pi*x[k])
    return y


def ref_shan(LB, UB, N, Fb, Fc):
    x = np.linspace(LB, UB, N)
    psi = np.sqrt(Fb)*(sinc2(Fb*x)*np.exp(2j*np.pi*Fc*x))
    return (psi, x)


def ref_fbsp(LB, UB, N, m, Fb, Fc):
    x = np.linspace(LB, UB, N)
    psi = np.sqrt(Fb)*((sinc2(Fb*x/m)**m)*np.exp(2j*np.pi*Fc*x))
    return (psi, x)


def ref_cmor(LB, UB, N, Fb, Fc):
    x = np.linspace(LB, UB, N)
    psi = ((np.pi*Fb)**(-0.5))*np.exp(2j*np.pi*Fc*x)*np.exp(-(x**2)/Fb)
    return (psi, x)


def ref_morl(LB, UB, N):
    x = np.linspace(LB, UB, N)
    psi = np.exp(-(x**2)/2)*np.cos(5*x)
    return (psi, x)


def ref_mexh(LB, UB, N):
    x = np.linspace(LB, UB, N)
    psi = (2/(np.sqrt(3)*np.pi**0.25))*np.exp(-(x**2)/2)*(1 - (x**2))
    return (psi, x)


def test_gaus():
    LB = -5
    UB = 5
    N = 1000
    for num in np.arange(1, 9):
        [psi, x] = ref_gaus(LB, UB, N, num)
        w = pywt.ContinuousWavelet("gaus" + str(num))
        PSI, X = w.wavefun(length=N)

        assert_allclose(np.real(PSI), np.real(psi))
        assert_allclose(np.imag(PSI), np.imag(psi))
        assert_allclose(X, x)


def test_cgau():
    LB = -5
    UB = 5
    N = 1000
    for num in np.arange(1, 9):
        [psi, x] = ref_cgau(LB, UB, N, num)
        w = pywt.ContinuousWavelet("cgau" + str(num))
        PSI, X = w.wavefun(length=N)

        assert_allclose(np.real(PSI), np.real(psi))
        assert_allclose(np.imag(PSI), np.imag(psi))
        assert_allclose(X, x)


def test_shan():
    LB = -20
    UB = 20
    N = 1000
    Fb = 1
    Fc = 1.5

    [psi, x] = ref_shan(LB, UB, N, Fb, Fc)
    w = pywt.ContinuousWavelet("shan{}-{}".format(Fb, Fc))
    assert_almost_equal(w.center_frequency, Fc)
    assert_almost_equal(w.bandwidth_frequency, Fb)
    w.upper_bound = UB
    w.lower_bound = LB
    PSI, X = w.wavefun(length=N)

    assert_allclose(np.real(PSI), np.real(psi), atol=1e-15)
    assert_allclose(np.imag(PSI), np.imag(psi), atol=1e-15)
    assert_allclose(X, x, atol=1e-15)

    LB = -20
    UB = 20
    N = 1000
    Fb = 1.5
    Fc = 1

    [psi, x] = ref_shan(LB, UB, N, Fb, Fc)
    w = pywt.ContinuousWavelet("shan{}-{}".format(Fb, Fc))
    assert_almost_equal(w.center_frequency, Fc)
    assert_almost_equal(w.bandwidth_frequency, Fb)
    w.upper_bound = UB
    w.lower_bound = LB
    PSI, X = w.wavefun(length=N)

    assert_allclose(np.real(PSI), np.real(psi), atol=1e-15)
    assert_allclose(np.imag(PSI), np.imag(psi), atol=1e-15)
    assert_allclose(X, x, atol=1e-15)


def test_cmor():
    LB = -20
    UB = 20
    N = 1000
    Fb = 1
    Fc = 1.5

    [psi, x] = ref_cmor(LB, UB, N, Fb, Fc)
    w = pywt.ContinuousWavelet("cmor{}-{}".format(Fb, Fc))
    assert_almost_equal(w.center_frequency, Fc)
    assert_almost_equal(w.bandwidth_frequency, Fb)
    w.upper_bound = UB
    w.lower_bound = LB
    PSI, X = w.wavefun(length=N)

    assert_allclose(np.real(PSI), np.real(psi), atol=1e-15)
    assert_allclose(np.imag(PSI), np.imag(psi), atol=1e-15)
    assert_allclose(X, x, atol=1e-15)

    LB = -20
    UB = 20
    N = 1000
    Fb = 1.5
    Fc = 1

    [psi, x] = ref_cmor(LB, UB, N, Fb, Fc)
    w = pywt.ContinuousWavelet("cmor{}-{}".format(Fb, Fc))
    assert_almost_equal(w.center_frequency, Fc)
    assert_almost_equal(w.bandwidth_frequency, Fb)
    w.upper_bound = UB
    w.lower_bound = LB
    PSI, X = w.wavefun(length=N)

    assert_allclose(np.real(PSI), np.real(psi), atol=1e-15)
    assert_allclose(np.imag(PSI), np.imag(psi), atol=1e-15)
    assert_allclose(X, x, atol=1e-15)


def test_fbsp():
    LB = -20
    UB = 20
    N = 1000
    M = 2
    Fb = 1
    Fc = 1.5

    [psi, x] = ref_fbsp(LB, UB, N, M, Fb, Fc)

    w = pywt.ContinuousWavelet("fbsp{}-{}-{}".format(M, Fb, Fc))
    assert_almost_equal(w.center_frequency, Fc)
    assert_almost_equal(w.bandwidth_frequency, Fb)
    w.fbsp_order = M
    w.upper_bound = UB
    w.lower_bound = LB
    PSI, X = w.wavefun(length=N)

    assert_allclose(np.real(PSI), np.real(psi), atol=1e-15)
    assert_allclose(np.imag(PSI), np.imag(psi), atol=1e-15)
    assert_allclose(X, x, atol=1e-15)

    LB = -20
    UB = 20
    N = 1000
    M = 2
    Fb = 1.5
    Fc = 1

    [psi, x] = ref_fbsp(LB, UB, N, M, Fb, Fc)
    w = pywt.ContinuousWavelet("fbsp{}-{}-{}".format(M, Fb, Fc))
    assert_almost_equal(w.center_frequency, Fc)
    assert_almost_equal(w.bandwidth_frequency, Fb)
    w.fbsp_order = M
    w.upper_bound = UB
    w.lower_bound = LB
    PSI, X = w.wavefun(length=N)

    assert_allclose(np.real(PSI), np.real(psi), atol=1e-15)
    assert_allclose(np.imag(PSI), np.imag(psi), atol=1e-15)
    assert_allclose(X, x, atol=1e-15)

    LB = -20
    UB = 20
    N = 1000
    M = 3
    Fb = 1.5
    Fc = 1.2

    [psi, x] = ref_fbsp(LB, UB, N, M, Fb, Fc)
    w = pywt.ContinuousWavelet("fbsp{}-{}-{}".format(M, Fb, Fc))
    assert_almost_equal(w.center_frequency, Fc)
    assert_almost_equal(w.bandwidth_frequency, Fb)
    w.fbsp_order = M
    w.upper_bound = UB
    w.lower_bound = LB
    PSI, X = w.wavefun(length=N)
    # TODO: investigate why atol = 1e-5 is necessary
    assert_allclose(np.real(PSI), np.real(psi), atol=1e-5)
    assert_allclose(np.imag(PSI), np.imag(psi), atol=1e-5)
    assert_allclose(X, x, atol=1e-15)


def test_morl():
    LB = -5
    UB = 5
    N = 1000

    [psi, x] = ref_morl(LB, UB, N)
    w = pywt.ContinuousWavelet("morl")
    w.upper_bound = UB
    w.lower_bound = LB
    PSI, X = w.wavefun(length=N)

    assert_allclose(np.real(PSI), np.real(psi))
    assert_allclose(np.imag(PSI), np.imag(psi))
    assert_allclose(X, x)


def test_mexh():
    LB = -5
    UB = 5
    N = 1000

    [psi, x] = ref_mexh(LB, UB, N)
    w = pywt.ContinuousWavelet("mexh")
    w.upper_bound = UB
    w.lower_bound = LB
    PSI, X = w.wavefun(length=N)

    assert_allclose(np.real(PSI), np.real(psi))
    assert_allclose(np.imag(PSI), np.imag(psi))
    assert_allclose(X, x)

    LB = -5
    UB = 5
    N = 1001

    [psi, x] = ref_mexh(LB, UB, N)
    w = pywt.ContinuousWavelet("mexh")
    w.upper_bound = UB
    w.lower_bound = LB
    PSI, X = w.wavefun(length=N)

    assert_allclose(np.real(PSI), np.real(psi))
    assert_allclose(np.imag(PSI), np.imag(psi))
    assert_allclose(X, x)


def test_cwt_parameters_in_names():

    for func in [pywt.ContinuousWavelet, pywt.DiscreteContinuousWavelet]:
        for name in ['fbsp', 'cmor', 'shan']:
            # additional parameters should be specified within the name
            assert_warns(FutureWarning, func, name)

        for name in ['cmor', 'shan']:
            # valid names
            func(name + '1.5-1.0')
            func(name + '1-4')

            # invalid names
            assert_raises(ValueError, func, name + '1.0')
            assert_raises(ValueError, func, name + 'B-C')
            assert_raises(ValueError, func, name + '1.0-1.0-1.0')

        # valid names
        func('fbsp1-1.5-1.0')
        func('fbsp1.0-1.5-1')
        func('fbsp2-5-1')

        # invalid name (non-integer order)
        assert_raises(ValueError, func, 'fbsp1.5-1-1')
        assert_raises(ValueError, func, 'fbspM-B-C')

        # invalid name (too few or too many params)
        assert_raises(ValueError, func, 'fbsp1.0')
        assert_raises(ValueError, func, 'fbsp1.0-0.4')
        assert_raises(ValueError, func, 'fbsp1-1-1-1')


def test_cwt_complex():
    for dtype in [np.float32, np.float64]:
        time, sst = pywt.data.nino()
        sst = np.asarray(sst, dtype=dtype)
        dt = time[1] - time[0]
        wavelet = 'cmor1.5-1.0'
        scales = np.arange(1, 32)

        # real-valued tranfsorm
        [cfs, f] = pywt.cwt(sst, scales, wavelet, dt)

        # complex-valued tranfsorm equals sum of the transforms of the real and
        # imaginary components
        [cfs_complex, f] = pywt.cwt(sst + 1j*sst, scales, wavelet, dt)
        assert_almost_equal(cfs + 1j*cfs, cfs_complex)


def test_cwt_small_scales():
    data = np.zeros(32)

    # A scale of 0.1 was chosen specifically to give a filter of length 2 for
    # mexh.  This corner case should not raise an error.
    cfs, f = pywt.cwt(data, scales=0.1, wavelet='mexh')
    assert_allclose(cfs, np.zeros_like(cfs))

    # extremely short scale factors raise a ValueError
    assert_raises(ValueError, pywt.cwt, data, scales=0.01, wavelet='mexh')

if __name__ == '__main__':
    run_module_suite()
