import numpy as np
from sklearn.mixture import BayesianGaussianMixture

X = np.array([[1, 2], [1, 4], [1, 0], [4, 2], [12, 4], [10, 7]])
print(X)
bgm = BayesianGaussianMixture(n_components=2, random_state=42).fit(X)
print(bgm)
bgm.means_
print(bgm.means_)

bgm.predict([[0, 0], [9, 3]])
print(bgm.predict([[0, 0], [9, 3]]))
