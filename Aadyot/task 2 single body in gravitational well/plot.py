import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("try3.txt");

plt.plot(0,0,'bo')
plt.plot(data["x"][0], data["y"][0], 'ro')
plt.plot(data["x"], data["y"], c='red')
plt.plot(data["x"].iloc[-1], data["y"].iloc[-1], 'ro')
plt.show()