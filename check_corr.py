import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Just run this with your arrays
def quick_check(arr1, arr2, arr3, arr4):
    sequences = [arr1, arr2, arr3, arr4]
    
    plt.figure(figsize=(10, 8))
    
    for i in range(4):
        plt.subplot(2, 2, i+1)
        
        for j in range(4):
            if i != j:
                corr = signal.correlate(sequences[i], sequences[j], mode='full')
                corr = corr / len(sequences[i])
                lags = np.arange(-len(sequences[i])+1, len(sequences[i]))
                plt.plot(lags, corr, alpha=0.7, label=f'vs Code {j+1}')
        
        plt.title(f'Code {i+1} cross-correlations')
        plt.xlabel('Lag')
        plt.ylabel('Correlation')
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

arr1 = np.array([0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1])
arr2 = np.array([0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0])
arr3 = np.array([0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0])
arr4 = np.array([0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1])

# Run it
quick_check(arr1, arr2, arr3, arr4)