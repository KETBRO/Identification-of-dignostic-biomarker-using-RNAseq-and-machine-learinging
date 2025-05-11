import pandas as pd
import numpy as np
from sklearn.model_selection import KFold, cross_val_score, GridSearchCV
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.linear_model import LassoCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from joblib import dump, load, parallel_backend

# Load normalized expression data and metadata
expression_data = pd.read_csv('normalized_expression_data.csv', index_col=0)
metadata = pd.read_csv('metadata.csv', index_col=0)

# Ensure metadata matches the sample order in the expression data
metadata = metadata.loc[expression_data.columns]


# Load DEGs from CSV
#degs_upregulated = pd.read_csv('upregulated_genes.csv')['x'].tolist()
#degs_downregulated = pd.read_csv('downregulated_genes.csv')['x'].tolist()

# Combine upregulated and downregulated genes
#degs = degs_upregulated + degs_downregulated

# Prepare data for machine learning
#X = expression_data.T[degs].values

# Prepare data for machine learning
X = expression_data.T.values  # Transpose to have samples as rows and genes as columns
y = metadata['description'].apply(lambda x: 1 if x == 'DR' else 0).values

# Scale the data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 5-Fold Cross-Validation Setup
kf = KFold(n_splits=5, shuffle=True, random_state=42)

# ------------------- SVM with RFE -------------------
svm = SVC(kernel="linear")
rfe = RFE(estimator=svm, n_features_to_select=21)
with parallel_backend('threading', n_jobs=-1):
    svm_scores = cross_val_score(rfe, X_scaled, y, cv=kf, n_jobs=-1)
rfe.fit(X_scaled, y)

# Correctly index the selected genes
#svm_selected_genes = [degs[i] for i in range(len(degs)) if rfe.support_[i]]
svm_selected_genes = expression_data.index[rfe.get_support()]

svm_selected_genes_df = pd.DataFrame(svm_selected_genes, columns=['svm_selected_genes'])
svm_selected_genes_df.to_csv('svm_selected_genes.csv', index=False)
# ------------------- LASSO Regression -------------------
lasso = LassoCV(cv=5, random_state=42, n_jobs=-1)
with parallel_backend('threading', n_jobs=-1):
    lasso_scores = cross_val_score(lasso, X_scaled, y, cv=kf, n_jobs=-1)
lasso.fit(X_scaled, y)

#lasso_selected_genes = [degs[i] for i in range(len(degs)) if lasso.coef_[i] != 0]
lasso_selected_genes = expression_data.index[lasso.coef_ != 0]

lasso_selected_genes_df = pd.DataFrame(lasso_selected_genes, columns=['lasso_selected_genes'])
lasso_selected_genes_df.to_csv('lasso_selected_genes.csv', index=False)

# ------------------- Random Forest -------------------
rf = RandomForestClassifier(random_state=42, n_jobs=-1)
param_grid_rf = {
    'n_estimators': [500, 1000, 2000],
    'max_features': [2, 5, 10, 'sqrt'],
    'max_depth': [None, 10, 20, 30]
}
rf_grid_search = GridSearchCV(rf, param_grid_rf, cv=kf, scoring='accuracy', n_jobs=-1)
with parallel_backend('threading', n_jobs=-1):
    rf_grid_search.fit(X_scaled, y)
best_rf_model = rf_grid_search.best_estimator_
rf_selected_genes_importances = best_rf_model.feature_importances_
indices = np.argsort(rf_selected_genes_importances)[-100:]  # Select top 30 genes based on importance
rf_selected_genes = expression_data.index[indices]

# Correctly index the selected genes
#rf_selected_genes = [degs[i] for i in indices]

rf_selected_genes_df = pd.DataFrame(rf_selected_genes, columns=['rf_selected_genes'])
rf_selected_genes_df.to_csv('rf_selected_genes_30.csv', index=False)


# ------------------- Gradient Boosting -------------------
gb = GradientBoostingClassifier(random_state=42)
param_grid_gb = {
    'n_estimators': [100, 200, 300],
    'learning_rate': [0.01, 0.1, 0.2],
    'max_depth': [3, 5, 7]
}
gb_grid_search = GridSearchCV(gb, param_grid_gb, cv=kf, scoring='accuracy', n_jobs=-1)
with parallel_backend('threading', n_jobs=-1):
    gb_grid_search.fit(X_scaled, y)
best_gb_model = gb_grid_search.best_estimator_
# Extract feature importances
feature_importances = best_gb_model.feature_importances_

# Get the indices of the most important features
top_indices = feature_importances.argsort()[::-1]  # Sort in descending order

# Get the top features (genes) - Adjust the number of top features as needed
top_n = 20  # For example, top 20 genes
top_genes = expression_data.index[top_indices[:top_n]]

# Save the top features (genes) to a CSV file
top_genes_df = pd.DataFrame(top_genes, columns=['Top Genes'])
top_genes_df.to_csv('top_gb_genes.csv', index=False)


common_genes = set(svm_selected_genes).intersection(set(lasso_selected_genes)).intersection(set(rf_selected_genes))#.intersection(set(top_genes))
common_genes_list = list(common_genes)

# Save common genes to a CSV file
common_genes_df = pd.DataFrame(common_genes_list, columns=['Common Genes'])
common_genes_df.to_csv('common_genes.csv', index=False)

#-------------------------------------------------------------------------------------------------

# # Print the best hyperparameters
# print("Best hyperparameters:", gb_grid_search.best_params_)

# # Evaluate the best model on the entire dataset (if needed)
# best_model_score = best_gb_model.score(X_scaled, y)
# print("Best model accuracy on the entire dataset:", best_model_score)

# # Save the best model to a file
# dump(best_gb_model, 'best_gb_model.joblib')

# # Save the best hyperparameters to a CSV file
# best_params_df = pd.DataFrame([gb_grid_search.best_params_])
# best_params_df.to_csv('best_gb_params.csv', index=False)

# print("Best model and hyperparameters saved to 'best_gb_model.joblib' and 'best_gb_params.csv'")

# # Extract feature importances
# feature_importances = best_gb_model.feature_importances_

# # Get the indices of the most important features
# top_indices = feature_importances.argsort()[::-1]  # Sort in descending order

# # Get the top features (genes) - Adjust the number of top features as needed
# top_n = 20  # For example, top 20 genes
# top_genes = expression_data.index[top_indices[:top_n]]

# # Save the top features (genes) to a CSV file
# top_genes_df = pd.DataFrame(top_genes, columns=['Top Genes'])
# top_genes_df.to_csv('top_gb_genes.csv', index=False)

# print("Top genes saved to 'top_gb_genes.csv'")


#-----------------------------------------------------------------------------------------------------




# ------------------- Identify common significant genes -------------------
common_genes = set(svm_selected_genes).intersection(set(lasso_selected_genes)).intersection(set(rf_selected_genes))
common_genes_list = list(common_genes)

# Save common genes to a CSV file
common_genes_df = pd.DataFrame(common_genes_list, columns=['Common Genes'])
common_genes_df.to_csv('common_genes.csv', index=False)

# Train final model using common significant genes
X_final = expression_data.T[common_genes_list].values
X_final_scaled = scaler.fit_transform(X_final)

# Ensemble Model: Combine Random Forest and Gradient Boosting
final_model_rf = RandomForestClassifier(n_estimators=2000, max_features=2, max_depth=20, random_state=42, n_jobs=-1)
final_model_rf.fit(X_final_scaled, y)

final_model_gb = GradientBoostingClassifier(n_estimators=300, learning_rate=0.1, max_depth=5, random_state=42)
final_model_gb.fit(X_final_scaled, y)

# Save the models and scaler
dump(final_model_rf, 'final_model_rf.joblib')
dump(final_model_gb, 'final_model_gb.joblib')
dump(scaler, 'scaler.joblib')

# Load the models and scaler for prediction
final_model_rf = load('final_model_rf.joblib')
final_model_gb = load('final_model_gb.joblib')
scaler = load('scaler.joblib')

# Load normalized expression data for a new sample
new_sample_data = pd.read_csv('new_sample_data.csv', index_col=0)
common_genes = pd.read_csv('common_genes.csv')['Common Genes']

# Ensure new sample data contains the common genes
new_sample_features = new_sample_data.T[common_genes].values

# Scale the new sample data
new_sample_scaled = scaler.transform(new_sample_features)

# Predict the class of the new sample using both models
prediction_rf = final_model_rf.predict(new_sample_scaled)
prediction_gb = final_model_gb.predict(new_sample_scaled)

# Ensemble Voting
final_prediction = np.round((prediction_rf + prediction_gb) / 2).astype(int)

# Output the prediction
if final_prediction[0] == 1:
    print("The sample is predicted to have DR.")
else:
    print("The sample is predicted to be non-DR.")
