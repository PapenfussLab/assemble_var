"""
Statistics module
"""

def mean(x):
    return float(sum(x))/len(x)


def median(x):
    y = copy.copy(x)
    y.sort()
    n = len(y)
    if n % 2: # Odd
        return y[n/2]
    else: # even
        return 0.5*(y[n/2]+y[n/2-1])
