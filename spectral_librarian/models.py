import numpy as np


class Polyfit1D(object):
    def __init__(self, x, y, w=None, deg=2):
        self.x = np.array(x)
        self.y = np.array(y)
        self.deg = deg
        self.w = None if w is None else np.array(w)

    def fit(self):
        self.coef = np.polyfit(x=self.x, y=self.y, deg=self.deg, w=self.w)

    def __call__(self, xx):
        return np.polyval(self.coef, xx)

    def fitted(self):
        return self(self.x)

    def res(self):
        return self.y - self(self.x)

    def plot(self, show=True):
        import matplotlib.pyplot as plt
        plt.plot(self.x, self.fitted(), c='red')
        if self.w is None:
            plt.scatter(self.x, self.y, c='orange')
        else:
            plt.scatter(self.x, self.y, c='orange', s=self.w / self.w.max() * 100)
        if show:
            plt.show()

    def error_l1(self):
        return np.abs(self.res()).sum()

    def average_l1(self):
        return np.abs(self.res()).sum()/len(self.x)


def polyfit(x, y, w=None, deg=1):
    """Fit a polynomial with least squares.

    Args:
        x (np.array): x coordinates
        y (np.array): y coordinates
        w (np.array): weights of observations
        deg (int): degree of the polynomial
    """
    P = Polyfit1D(x, y, w, deg)
    P.fit()
    return P