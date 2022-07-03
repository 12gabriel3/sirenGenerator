from scipy import signal
from scipy.io.wavfile import write
from numpy import arange, float64, pi, int16, append, exp, cumsum, flip
import matplotlib.pyplot as plt


def discharge(Vs, t, RC):
    return (Vs) * exp(-t/(RC)) 

def charge(Vs, V0, t, RC):
    return (Vs - V0) * (1 - exp(-t/(RC))) + V0

def cycle(N, V0, Tdischarge, Tcharge, ExtraCycles):
    RCdis = Tdischarge / 2
    RCcha = Tcharge / 2
    chargeEnd = charge(2, V0, Tcharge , RCcha)
    dischargeEnd = discharge(chargeEnd, Tdischarge, RCdis)
    if(ExtraCycles):
        return append(charge(2, V0, arange(Tcharge * N) / N, RCcha), \
            append(discharge(chargeEnd, arange(Tdischarge * N) / N, RCdis), \
            cycle(N, dischargeEnd, Tdischarge, Tcharge, ExtraCycles - 1)))
    else:
        return append(charge(2, V0, arange(Tcharge * N) / N, RCcha), discharge(chargeEnd, arange(Tdischarge * N) / N, RCdis))

def f(f_c, N, V0, Tcharge, Tdischarge, Ncycles):
    if(not Ncycles):
        return 0
    sum = cumsum(cycle(N, V0, Tdischarge, Tcharge, Ncycles - 1).astype(float64) - 1)
    t = arange(N * (Tdischarge + Tcharge) * Ncycles) / N
    if t.size > sum.size:
        t = t[:sum.size]
    elif t.size < sum.size:
        sum = sum[:t.size]
    return signal.sawtooth(2*pi*f_c*t - (sum * 1 / N * -5000), 0.5)

def to_integer(signal):
    # Take samples in [-1, 1] and scale to 16-bit integers,
    # values between -2^15 and 2^15 - 1.
    return int16(signal*(2**15 - 1))

def ClampZeroCrossing(signal):
    i = next(x[0] for x in enumerate(signal) if abs(x[1]) < 100)
    flipped = flip(signal)

    if signal[i] > signal[i - 1]:
        j = next(i for i, v in enumerate(flipped) if abs(v) < 100 and v < flipped[i - 1])
    else:
        j = next(i for i, v in enumerate(flipped) if abs(v) < 100 and v > flipped[i - 1])

    return signal[i:signal.size - j]

N = 44100 # samples per secon

data = f(1100, N, 0, 0.1, 1, 1)
write("siren.wav", N, ClampZeroCrossing(to_integer(data)))