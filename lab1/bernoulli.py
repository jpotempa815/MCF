import numpy as np
import matplotlib.pyplot as plt

def load_data(filename):
    data = {0.1: [], 0.5: [], 0.9: []}
    p_values = [0.1, 0.5, 0.9]
    
    with open(filename, 'r') as file:
        current_p = 0.1
        for line in file:
            values = list(map(float, line.split()))
            data[current_p].append(values)
            
            if len(data[current_p]) == 6:  
                index = p_values.index(current_p)
                if index < len(p_values) - 1:
                    current_p = p_values[index + 1]
    
    for key in data:
        data[key] = np.array(data[key])
    return data

def plot_errors(data):
    plt.figure(figsize=(12, 5))
    
    for i, (p, values) in enumerate(data.items()):
        n, _, err_X, err_var = values.T
        
        plt.subplot(1, 2, 1)
        plt.loglog(n, err_X, marker='o', label=f'p={p}')
        plt.xlabel('Liczba losowań N')
        plt.ylabel('Błąd względny wartości oczekiwanej')
        plt.legend()
        plt.grid(True, which='both', linestyle='--')
        
        plt.subplot(1, 2, 2)
        plt.loglog(n, err_var, marker='s', label=f'p={p}')
        plt.xlabel('Liczba losowań N')
        plt.ylabel('Błąd względny wariancji')
        plt.legend()
        plt.grid(True, which='both', linestyle='--')
    
    plt.tight_layout()
    plt.savefig("bernoulli_errors.png", dpi=300)

data = load_data("lab1.txt")
plot_errors(data)
