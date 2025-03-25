import matplotlib.pyplot as plt 
import numpy as np 
from scipy.stats import chi2

k = 10

def dystrybuanta(x_k, x_0):
    prop = 4/5 * ((x_k + x_k**2/2 - x_k**4/4) - (x_0 + x_0**2/2 - x_0**4/4))
    return prop

fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(nrows=2, ncols=2)

def histogram(ax, nazwa, nr):
    data = np.loadtxt(f"{nazwa}.txt")
    counts, bins, _ = ax.hist(data, bins=k, alpha=0.6, density=True, color='plum', edgecolor='black')
    x = np.linspace(0,1, 100)
    fx = (4/5) * (1 + x - x**3)
    ax.plot(x, fx, lw=2, label='fgp', color='mediumvioletred')
    ax.set_title(f"{nr})")
    
    chi_square_test(data, counts, bins, nazwa)
    # print(chi_kwadrat(data, counts, bins))

def chi_square_test(data, counts, bins, nazwa):
    N = len(data)
    expected = np.array([dystrybuanta(bins[i+1], bins[i]) * N for i in range(k)])
    observed = counts * np.diff(bins) * N
    chi_stat = np.sum((observed - expected) ** 2 / (expected))
    
    critical_value = chi2.ppf(0.95, df=k-1)
    print(f"Test chi-kwadrat dla {nazwa}: statystyka = {chi_stat:.3f}, wartość krytyczna = {critical_value:.3f}")
    if chi_stat < critical_value:
        print(f"Hipoteza zerowa dla {nazwa} nie została odrzucona. Dane pasują do rozkładu.")
    else:
        print(f"Hipoteza zerowa dla {nazwa} została odrzucona. Dane mogą nie pasować do rozkładu.")

def chi_kwadrat(data, counts, bins):
    N = len(data)
    chi = 0
    for i in range(k):
        p_i = dystrybuanta(bins[i+1], bins[i])
        n_i = counts[i] * 0.1 * N
        chi += (n_i - p_i*N)**2 / (p_i*N)
    return chi

histogram(ax0, "zlozony", 1)
histogram(ax1, "lancuch_1", 2)
histogram(ax2, "lancuch_2", 3)
histogram(ax3, "eliminacja", 4)

fig.tight_layout()
plt.savefig("histogram.png")