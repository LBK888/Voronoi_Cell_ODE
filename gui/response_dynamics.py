"""Rise tangent and peak-interval analysis of uniformly sampled spatial means."""

import numpy as np
from scipy.signal import find_peaks


def response_dynamics(values, dt):
    """Return NaN for features the sampled trace cannot resolve."""
    slope = onset = frequency = np.nan
    y = np.asarray(values, dtype=float)
    if len(y) >= 5 and np.isfinite(y).all():
        # A five-frame local linear regression estimates the tangent without
        # differentiating adjacent noisy samples. Its fitted midpoint is the mean.
        centers = np.convolve(y, np.ones(5) / 5, mode='valid')
        slopes = np.correlate(y - y[0], np.arange(-2., 3.), mode='valid') / (10 * dt)
        smooth = np.r_[y[0], centers, y[-1]]
        noise = np.median(np.abs(np.diff(y, n=2))) / 1.652
        prominence = max(0.10 * np.ptp(smooth), 5 * noise,
                         64 * np.finfo(float).eps * max(np.max(np.abs(y)), 1.))
        peaks = find_peaks(smooth, prominence=prominence)[0]
        troughs = find_peaks(-smooth, prominence=prominence)[0]
        # Endpoint inclusion also admits a monotone or unfinished rise.
        for end in np.r_[peaks, len(smooth) - 1]:
            previous = troughs[troughs < end]
            start = previous[-1] if previous.size else 0
            if smooth[end] - smooth[start] <= prominence:
                continue
            # Interior smooth index i corresponds to raw sample i+1.
            lo = start if start else 0
            hi = min(end, len(slopes))
            if hi <= lo:
                continue
            index = lo + np.argmax(slopes[lo:hi])
            if slopes[index] > 0:
                slope = float(slopes[index])
                start_time = (start + 1) * dt if start else 0.
                tangent = (index + 2) * dt - (centers[index] - smooth[start]) / slope
                onset = float(np.clip(tangent, start_time, (index + 2) * dt))
                break
        if len(peaks) >= 3:
            intervals = np.diff(peaks) * dt
            period = np.median(intervals)
            if np.all(np.abs(intervals - period) <= 0.25 * period) and period >= 5 * dt:
                lag = round(period / dt)
                left, right = y[:-lag], y[lag:]
                if np.ptp(left) and np.ptp(right) and np.corrcoef(left, right)[0, 1] >= 0.7:
                    frequency = float(1 / period)
    return dict(rise_slope=slope, onset_time=onset, frequency=frequency)
