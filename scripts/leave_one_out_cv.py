


import sklearn
sklearn.model_selection.cross_val_score(model, X, y, scoring = 'r2')


import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split, cross_val_score, LeaveOneOut
from sklearn import metrics

admit = pd.read_csv('assets/admissions.csv')
admit = admit.dropna()
Xr = admit[['admit', 'gre', 'prestige']]
yr = admit['gpa']
X_array = np.array(Xr) #r stands for 'regression'
y_array = np.array(yr)

scores = cross_val_score(LinearRegression(), Xr, yr, cv=397, scoring = "r2")
print("Cross-validated scores:", scores)
print("Average: ", scores.mean())
print("Variance: ", np.std(scores))
