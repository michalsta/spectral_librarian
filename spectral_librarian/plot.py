def plot_spectrum(mz, i, show=True):
    import matplotlib.pyplot as plt
    plt.stem(mz, i)
    if show:
        plt.show()
