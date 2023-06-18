import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("try1.txt");
plt.plot(0,0,'bo')
plt.plot(data["x"][0], data["y"][0], 'ro')
plt.plot(data["x"], data["y"], c='red')
plt.show()