def plot_spectrum(mz, i, show=True):
    import matplotlib.pyplot as plt
    plt.stem(mz, i)
    if show:
        plt.show()


def side_plot_spectra(spec0, spec1):
    import matplotlib.pyplot as plt
    plt.stem(spec0[0], spec0[1])
    plt.stem(spec1[0],-spec1[1])
    if show:
        plt.show()

