import matplotlib.pyplot as plt
import numpy as np

# Read the numbers from the file
with open('distrib.txt', 'r') as file:
    numbers = [float(line.strip()) for line in file]

# Plot the distribution
plt.hist(numbers, bins=np.logspace(np.log10(min(numbers)), np.log10(max(numbers)), 50))
plt.xscale('log')

# Set plot title and labels
plt.title('Number Distribution')
plt.xlabel('Numbers (Log Scale)')
plt.ylabel('Frequency')

# Show the plot
plt.show()
