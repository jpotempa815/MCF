import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def wczytywanie(file):
    data = []
    for i in file:
        data = np.loadtxt(file)
    return data

def rys(nazwa, lit):
    data1 = wczytywanie(f"{nazwa}.txt")
    df = pd.DataFrame(data1)
    plt.subplots(figsize=(5, 5))
    plt.grid(linestyle='--')
    plt.scatter(df[0], df[1], s=2, marker='.', color='m')
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(f"{lit})")
    plt.savefig(f"{nazwa}.png")

rys('box_muller', 'a')
rys('kolo', 'b')
rys('elip', 'b')
rys('elip_norm', 'a')