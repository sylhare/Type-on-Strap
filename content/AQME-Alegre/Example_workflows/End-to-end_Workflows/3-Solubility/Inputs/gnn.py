import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from gnn_functions import *
import sklearn.metrics as metrics
import tensorflow as tf


sol = pd.read_csv('solubility_xtb.csv')
valid, test, train = np.split(sol[['smiles','xtbjson']].sample(frac=1., random_state=100), [50, 100])


train_dataset, valid_dataset, test_dataset = gnn_data(valid, test, train, sol)
inputs, outputs = next(train_dataset.as_numpy_iterator())


model = gnn_model()
model.compile(loss='mae', optimizer=tf.keras.optimizers.Adam(1E-3))
print('training')
model.fit(train_dataset, validation_data=valid_dataset, epochs=500)


# Predict solubility of the external test set
test_predictions = model.predict(test_dataset)
test_db_values = sol.set_index('smiles').reindex(test.smiles)['measured log solubility in mols per litre'].values

print('testing')
# Plot the results
fig = plt.subplots(figsize=(3,3))

ax1 = sns.scatterplot(test_db_values,test_predictions.flatten(),s=30,marker='o',color='b',alpha=0.5)
ax1.set_xlabel(r'Measured',fontsize=10)
ax1.set_ylabel(r'Predicted',fontsize=10)
# ax1.grid(linestyle='--', linewidth=1)
mae = metrics.mean_absolute_error(test_db_values,test_predictions.flatten())
r2 = metrics.r2_score(test_db_values,test_predictions.flatten())

plt.annotate(f"$R^2$ = {round(r2,1)} \nMAE = {round(mae,1)} ", xy=(-2.5, -5.6), fontsize=10)
plt.savefig('solubility-gnn.jpg',dpi=400)
plt.show()




