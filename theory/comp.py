import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('/home/kupalnic/4th_sem/System_of_linear_equations/build/benchmark_results.csv')

for s in df['Sparsity'].unique():
    subset = df[df['Sparsity'] == s]
    plt.figure(figsize=(10, 6))
    plt.plot(subset['Size'], subset['Dense_ms'], label='Matrix', marker='o')
    plt.plot(subset['Size'], subset['CSR_ms'], label='CSR Matrix', marker='s')
    
    plt.title(f'Сравнение скорости при разреженности {s*100}%')
    plt.xlabel('Размер матрицы (N x N)')
    plt.ylabel('Время (микросекунды)')
    plt.legend()
    plt.grid(True)
    plt.yscale('log') 
    plt.show()
