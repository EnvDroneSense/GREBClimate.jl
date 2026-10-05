# Recorded results that test/test_golden.jl compares a run on the NCEP dataset
# with. Monthly unweighted global means of Ts and Ta (K) and q (kg/kg).
# A model change made on purpose means replacing these values.

# :full_model on the stored flux corrections, 1 control year
const GOLDEN_CTRL = [
    (Ts = 276.6376, Ta = 279.01324, q = 0.006483287),
    (Ts = 276.08505, Ta = 278.59598, q = 0.0066876207),
    (Ts = 275.80753, Ta = 278.15857, q = 0.0068251393),
    (Ts = 276.7514, Ta = 278.90686, q = 0.0069842925),
    (Ts = 278.5325, Ta = 280.68484, q = 0.007306205),
    (Ts = 280.12466, Ta = 282.42026, q = 0.007825737),
    (Ts = 280.6912, Ta = 283.14557, q = 0.008249164),
    (Ts = 280.36264, Ta = 282.90518, q = 0.008242705),
    (Ts = 279.22748, Ta = 281.76135, q = 0.007820126),
    (Ts = 278.26984, Ta = 280.7673, q = 0.0074224365),
    (Ts = 277.96216, Ta = 280.53217, q = 0.0072685555),
    (Ts = 277.9465, Ta = 280.61523, q = 0.007346397),
]

# The same run, 1 scenario year (anomaly against the control)
const GOLDEN_SCNR = [
    (Ts = -0.0076227454, Ta = -0.0071252817, q = -4.826931e-08),
    (Ts = -0.0013025337, Ta = -0.0015087194, q = 1.1220921e-08),
    (Ts = -0.00011379851, Ta = -0.0001452234, q = -7.757092e-09),
    (Ts = 0.000121321944, Ta = 0.000109407636, q = 1.9124322e-09),
    (Ts = 0.00016302532, Ta = 0.00015985966, q = 1.4248567e-08),
    (Ts = 9.6678734e-05, Ta = 0.000101053054, q = 1.7348425e-08),
    (Ts = 5.298853e-05, Ta = 5.4988595e-05, q = 1.2192554e-08),
    (Ts = 3.7090645e-05, Ta = 3.6418438e-05, q = 5.587703e-09),
    (Ts = 8.727444e-05, Ta = 8.031726e-05, q = 5.9890044e-09),
    (Ts = 0.00010110272, Ta = 0.00010247363, q = 2.6633937e-09),
    (Ts = 7.9893405e-05, Ta = 8.220805e-05, q = -1.2187116e-09),
    (Ts = 5.7422454e-05, Ta = 5.850858e-05, q = -2.5535258e-09),
]

# :full_model, 1 control year after a one-year spin-up
const GOLDEN_FLUX = [
    (Ts = 277.28638, Ta = 279.7823, q = 0.0063898247),
    (Ts = 276.35574, Ta = 278.86285, q = 0.0064106397),
    (Ts = 275.92465, Ta = 278.2391, q = 0.006436156),
    (Ts = 276.80832, Ta = 278.9109, q = 0.0065385858),
    (Ts = 278.5104, Ta = 280.59595, q = 0.006805867),
    (Ts = 279.90192, Ta = 282.09723, q = 0.0072038374),
    (Ts = 280.28174, Ta = 282.58386, q = 0.007454401),
    (Ts = 279.82877, Ta = 282.20816, q = 0.007350458),
    (Ts = 278.68015, Ta = 281.04367, q = 0.0068896553),
    (Ts = 277.73822, Ta = 280.0452, q = 0.0064565144),
    (Ts = 277.4416, Ta = 279.80316, q = 0.006277522),
    (Ts = 277.36218, Ta = 279.8118, q = 0.0062873242),
]

# MSCM (Monash Simple Climate Model) database, 2xCO2 with every process on:
# year-1 area-weighted global-mean surface temperature response, K
const MSCM_YEAR1 = 0.594636
