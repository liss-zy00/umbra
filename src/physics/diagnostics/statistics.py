from __future__ import annotations

"""
Light-weight correlation / radial averaging helpers for pipeline metrics.

Used to summarize outer-disk kinematics where non-baryonic effects dominate in SPARC profiles.
"""

import numpy as np


def edge_radial_mean(values, fraction: float = 0.2, min_points: int = 3) -> float | None:
    """
    Calculates the mean of the outer portion (edge) of a galactic disk.
    
    Logic:
    - Dark matter (or our magnetic pressure equivalent) dominates the outer edges of galaxies where 
      baryonic matter drops off. Therefore, we evaluate the gap specifically at the outer 'fraction'.
    - We extract the last `fraction` (e.g., 20%) of the valid data points.
    - We enforce a minimum number of points (min_points) to avoid statistical noise from a single stray pixel.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    
    if arr.size == 0:
        return None
        
    # Determine how many points constitute the 'edge'. 
    # It takes the max between the defined minimum points and the calculated fraction.
    n = min(arr.size, max(min_points, int(np.ceil(arr.size * fraction))))
    
    # Return the mean of the last 'n' points (the outermost radial measurements)
    return float(np.mean(arr[-n:]))


def pearson_correlation_coefficient(a: np.ndarray, b: np.ndarray) -> float | None:
    """
    Calculates the Pearson correlation coefficient between two arrays.
    
    Logic:
    - This is used to check if the curve of the Magnetic Acceleration follows the exact 
      shape/trend of the Acceleration Gap across the radius. A value close to 1.0 means perfect correlation.
    - It filters out NaN or infinite values first.
    - Requires at least 3 valid points to make a mathematically sound correlation.
    - If the standard deviation is zero (a flat horizontal line), correlation is mathematically undefined, returning None.
    """
    valid = np.isfinite(a) & np.isfinite(b)
    if np.count_nonzero(valid) < 3:
        return None
        
    a_v = a[valid]
    b_v = b[valid]
    
    if np.std(a_v) == 0 or np.std(b_v) == 0:
        return None
        
    # Return the off-diagonal correlation value from the correlation matrix
    return float(np.corrcoef(a_v, b_v)[0, 1])