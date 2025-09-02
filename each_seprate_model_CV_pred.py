# To ensure robustness and to avoid information leakage, feature selection was performed only on the training set. 
# Several feature selection strategies were employed independently, including Support Vector Machine Recursive Feature Elimination (SVM-RFE), 
# Least Absolute Shrinkage and Selection Operator (LASSO), Mann–Whitney U test, Boruta algorithm, distance correlation (DCor), random forest 
# feature importance, and Elastic Net regularization. For each method, the number of top-ranked features was progressively reduced from ten to five, 
# in steps of one, and classifiers were retrained at each step to evaluate performance. Among all methods, the Elastic Net–based approach identified 
# a seven-gene panel that yielded the best diagnostic accuracy.
A=2
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold, cross_val_score, GridSearchCV, train_test_split, cross_val_predict
from sklearn.svm import SVC
from sklearn.feature_selection import RFE, SelectFromModel
from sklearn.linear_model import LassoCV, LogisticRegression, LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import make_scorer, mean_squared_error, accuracy_score, matthews_corrcoef, roc_auc_score, precision_score, recall_score, f1_score, confusion_matrix, roc_curve, auc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from joblib import dump, load, parallel_backend
from sklearn.pipeline import Pipeline
#import warnings
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from scipy.stats import mannwhitneyu
from boruta import BorutaPy
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata
from sklearn.linear_model import ElasticNetCV
import joblib


# Load normalized expression data and metadata
expression_data = pd.read_csv('normalized_expression_data.csv', index_col=0)
metadata = pd.read_csv('metadata.csv', index_col=0)

# Ensure metadata matches the sample order in the expression data
metadata = metadata.loc[expression_data.columns]

# Load DEGs from CSV
degs_upregulated = pd.read_csv('upregulated_genes.csv')['x'].tolist()
degs_downregulated = pd.read_csv('downregulated_genes.csv')['x'].tolist()

# Combine upregulated and downregulated genes
degs = degs_upregulated + degs_downregulated

Z = expression_data.T
X = Z[degs] 

# Prepare data for machine learning
#X = expression_data.T[degs].values
y = metadata['description'].apply(lambda x: 1 if x == 'DR' else 0).values



# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

# Scale the data
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)
#--------------------------------------------------------------------------------------------------------------------
# # Define the model
# svm = SVC(kernel="linear", probability=True)

# # Setup RFE with SVM, initially without specifying the number of features
# rfe = RFE(estimator=svm)

# # Create a pipeline
# pipeline = Pipeline([
#     ('feature_selection', rfe),
#     ('classification', svm)
# ])

# # Define a parameter grid
# param_grid = {
#     'feature_selection__n_features_to_select': [8, 9, 10],  # Testing different numbers of features
#     'classification__C': [0.1, 1, 10, 100]  # Regularization parameter
# }

# # Create a scoring function
# scoring = make_scorer(accuracy_score)
kf = KFold(n_splits=5, shuffle=True, random_state=42)
# # Setup the GridSearchCV
# grid_search = GridSearchCV(pipeline, param_grid, cv=kf, scoring=scoring, verbose=1)

# # Fit GridSearchCV
# SVM_rfe=grid_search.fit(X_train_scaled, y_train)
# best_model = SVM_rfe.best_estimator_
# # score=cross_val_score(best_model, X_train_scaled, y_train, cv=kf)
# # print("cross_val_score", score.mean())
#--------------------------------------------------------------
# test_score = best_model.score(X_test_scaled, y_test)
# print("Test score: ", test_score)
#----------------------------------------------------------------
# selected_features_mask = best_model.named_steps['feature_selection'].support_
# selected_features = low_corr_data.columns[selected_features_mask]

# X_train_scaled_selected_svm_rfe = X_train_scaled[:, selected_features_mask]
# X_test_scaled_selected_svm_rfe = X_test_scaled[:, selected_features_mask]

# X_train_selected_svm_rfe = X_train[:, selected_features_mask]
# X_test_selected_svm_rfe = X_test[:, selected_features_mask]

# selected_indices = np.where(selected_features_mask)[0]
# X_train_selected_svm_rfe = X_train[:, selected_indices]
# X_test_selected_svm_rfe = X_test[:, selected_indices]

# Assuming the best SVM settings are stored within the pipeline
# svm_classifier = best_model.named_steps['classification']
# svm_classifier.fit(X_train_selected_svm_rfe, y_train)
# score=cross_val_score(best_model, X_train_selected_svm_rfe, y_train, cv=kf)
# print("cross_val_score", score.mean())
#----------------------------------------------------------------------------------
# # Perform LassoCV for feature selection
# lasso = LassoCV(alphas=np.logspace(-3, 3, 10), cv=5, max_iter=10000)
# lasso.fit(X_train_scaled, y_train)

# # Sort features by their absolute coefficients
# lasso_coef = np.abs(lasso.coef_)
# sorted_features = np.argsort(lasso_coef)[::-1]

# # Select the top N features
# N = 10  # Number of top features to select(from all feature(498) give better result)
# top_features = sorted_features[:N]

# # Transform data to include only the top N features
# X_train_scaled_lasso_selected = X_train_scaled[:, top_features]
# X_test_scaled_lasso_selected = X_test_scaled[:, top_features]

# #-----------------------------------------------------------------------------------

# # Perform Mann-Whitney U test for each feature
# p_values = []
# for i in range(X_train_scaled.shape[1]):
#     # Perform Mann-Whitney U test between the two classes for each feature
#     u_stat, p_val = mannwhitneyu(X_train_scaled[y_train == 0, i], X_train_scaled[y_train == 1, i])
#     p_values.append(p_val)

# # Convert p-values to a numpy array and sort features by p-value
# p_values = np.array(p_values)
# sorted_features = np.argsort(p_values)

# # Select the top N features
# N = 10  # Number of top features to select
# top_features = sorted_features[:N]

# # Transform data to include only the top N features
# X_train_scaled_mannwhitney_selected = X_train_scaled[:, top_features]
# X_test_scaled_mannwhitney_selected = X_test_scaled[:, top_features]
# #---------------------------------------------------------------------------------

# # Define the random forest classifier
# rf = RandomForestClassifier(n_estimators=100, random_state=42)

# # Define Boruta feature selection method
# boruta_selector = BorutaPy(rf, n_estimators='auto', random_state=42)

# # Fit the Boruta selector
# boruta_selector.fit(X_train_scaled, y_train)

# # Get the selected features mask
# selected_features_mask = boruta_selector.support_

# # Get the ranking of the features
# feature_ranking = boruta_selector.ranking_
# sorted_features = np.argsort(feature_ranking)
# n_top_features = 10
# top_features = sorted_features[:n_top_features]


# # Transform data to include only the top N features
# X_train_scaled_boruta_selected = X_train_scaled[:, top_features]
# X_test_scaled_boruta_selected = X_test_scaled[:, top_features]

#-------------------------------------------------------------------------------------

# # Define a function to calculate distance correlation
# def distance_correlation(X, y):
#     def dist_corr(X, Y):
#         """ Compute the distance correlation function """
#         n = X.shape[0]
#         a = squareform(pdist(X.reshape(-1, 1)))
#         b = squareform(pdist(Y.reshape(-1, 1)))

#         A = a - a.mean(axis=0)[None, :] - a.mean(axis=1)[:, None] + a.mean()
#         B = b - b.mean(axis=0)[None, :] - b.mean(axis=1)[:, None] + b.mean()

#         dCovXY = np.sqrt((A * B).sum() / (n ** 2))
#         dVarXX = np.sqrt((A * A).sum() / (n ** 2))
#         dVarYY = np.sqrt((B * B).sum() / (n ** 2))

#         return dCovXY / np.sqrt(dVarXX * dVarYY)

#     return np.array([dist_corr(X[:, i], y) for i in range(X.shape[1])])

# # Perform DCor for feature selection
# dcor_values = distance_correlation(X_train_scaled, y_train)

# # Sort features by their DCor values
# sorted_features = np.argsort(dcor_values)[::-1]

# # Select the top N features
# N = 9  # Number of top features to select
# top_features = sorted_features[:N]

# # Transform data to include only the top N features
# X_train_scaled_dcor_selected = X_train_scaled[:, top_features]
# X_test_scaled_dcor_selected = X_test_scaled[:, top_features]

# #------------------------------------------------------------------------------
# # Train a Random Forest model
# rf = RandomForestClassifier(n_estimators=100, random_state=42)
# rf.fit(X_train_scaled, y_train)

# # Extract feature importances
# importances = rf.feature_importances_

# # Sort features by their importance
# sorted_indices = np.argsort(importances)[::-1]

# # Select the top N features
# N = 5  # Number of top features to select
# top_features = sorted_indices[:N]

# # Transform data to include only the top N features
# X_train_scaled_rf_selected = X_train_scaled[:, top_features]
# X_test_scaled_rf_selected = X_test_scaled[:, top_features]

#------------------------------------------------------------------------------

# Train ElasticNetCV for feature selection
elastic_net = ElasticNetCV(l1_ratio=[.1, .5, .7, .9, .95, .99, 1], alphas=np.logspace(-1, 1, 10), cv=5, max_iter=20000)
elastic_net.fit(X_train_scaled, y_train)

# Sort features by their absolute coefficients
elastic_net_coef = np.abs(elastic_net.coef_)
sorted_features = np.argsort(elastic_net_coef)[::-1]

# Select the top N features
N = 7  # Number of top features to select
top_features = sorted_features[:N]

# Get feature names from original DataFrame
feature_names = X_train.columns
selected_feature_names = feature_names[top_features]

# Save selected features
selected_features_df = pd.DataFrame({
    'Feature': selected_feature_names,
    'Coefficient': elastic_net_coef[top_features]
})


#selected_features_df.to_csv('selected_features_elastic_net_6.csv', index=False)

# Transform data to include only the top N features
X_train_scaled_enet_selected = X_train_scaled[:, top_features]
X_test_scaled_enet_selected = X_test_scaled[:, top_features]

#------------------------------------------------------------------------------
# Define models
models = {
    'Logistic Regression': LogisticRegression(),
    'Random Forest': RandomForestClassifier(n_estimators=100),
    'Gradient Boosting': GradientBoostingClassifier(n_estimators=100),
    'SVM': SVC(probability=True),
    'Neural Network': MLPClassifier(hidden_layer_sizes=(100,), max_iter=10000)
}

# Evaluate each model using cross-validation
for name, model in models.items():
    scores = cross_val_score(model, X_train_scaled_enet_selected, y_train, cv=kf, scoring='accuracy')
    print(f"{name} CV_Accuracy: {np.mean(scores):.4f} ± {np.std(scores):.4f}")
    
   # Calculate AUC for each fold
    y_proba = cross_val_predict(model, X_train_scaled_enet_selected, y_train, cv=kf, method='predict_proba')
    auc_scores = []
    for train_idx, test_idx in kf.split(X_train_scaled_enet_selected):
        fold_auc = roc_auc_score(y_train[test_idx], y_proba[test_idx, 1])
        auc_scores.append(fold_auc)
        print(f"{name} Fold AUC: {fold_auc:.4f}")
    
    mean_auc = np.mean(auc_scores)
    print(f"{name} Mean AUC: {mean_auc:.4f}")    
#-----------------------------------------------------------------------------------    
#     # Train and evaluate models on the test set
#     results = {}
#     for name, model in models.items():
#         # Fit model
#         model.fit(X_train_scaled_enet_selected, y_train)
        
#         # Predict and evaluate
#         y_pred = model.predict(X_test_scaled_enet_selected)
#         y_proba_ev = model.predict_proba(X_test_scaled_enet_selected)[:, 1]
        
#         accuracy = accuracy_score(y_test, y_pred)
#         auc = roc_auc_score(y_test, y_proba_ev)
#         results[name] = (accuracy, auc)

#     # Print results
#     for name, (acc, auc) in results.items():
#         print(f"{name} - Accuracy: {acc:.2f}, AUC: {auc:.2f}")
        
    
# #-----------------------------------------------------------------------------

# Train and evaluate models on the test set
results = {}
for name, model in models.items():
    # Fit model
    model.fit(X_train_scaled_enet_selected, y_train)
    
    # Predict and evaluate
    y_pred = model.predict(X_test_scaled_enet_selected)
    y_proba_ev = model.predict_proba(X_test_scaled_enet_selected)[:, 1]
    
    accuracy = accuracy_score(y_test, y_pred)
    mse = mean_squared_error(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    cm = confusion_matrix(y_test, y_pred)
    tn, fp, fn, tp = cm.ravel()
    sensitivity = recall  # Same as recall
    specificity = tn / (tn + fp)
    mcc = matthews_corrcoef(y_test, y_pred)
    auc = roc_auc_score(y_test, y_proba_ev)
    
    results[name] = {
        'Accuracy': accuracy,
        'MSE': mse,
        'F1 Score': f1,
        'Precision': precision,
        'Recall': recall,
        'Sensitivity': sensitivity,
        'Specificity': specificity,
        'MCC': mcc,
        'AUC': auc
    }

# Print results
for name, metrics in results.items():
    print(f"{name}:")
    for metric, value in metrics.items():
        print(f"  {metric}: {value:.4f}")



#-----------------------------------------------------------------------------

# Select SVM model 
svm_model = models['SVM'] #(Chnage the model seprate figure)

# Calculate AUC for each fold and store FPR and TPR for each fold
fprs = []
tprs = []
auc_scores = []

for train_idx, test_idx in kf.split(X_train_scaled_enet_selected):
    svm_model.fit(X_train_scaled_enet_selected[train_idx], y_train[train_idx])
    y_proba = svm_model.predict_proba(X_train_scaled_enet_selected[test_idx])[:, 1]
    fpr, tpr, _ = roc_curve(y_train[test_idx], y_proba)
    fprs.append(fpr)
    tprs.append(tpr)
    fold_auc = roc_auc_score(y_train[test_idx], y_proba)
    auc_scores.append(fold_auc)
    print(f"SVM Fold AUC: {fold_auc:.4f}")

mean_auc = np.mean(auc_scores)
print(f"SVM Mean AUC: {mean_auc:.4f}")

# Plot ROC curves for each fold
plt.figure()
for i in range(len(fprs)):
    plt.plot(fprs[i], tprs[i], lw=2, label=f'Fold {i+1} ROC curve (AUC = {auc_scores[i]:.2f})')

# Plot mean ROC curve
mean_fpr = np.linspace(0, 1, 100)
mean_tpr = np.mean([np.interp(mean_fpr, fprs[i], tprs[i]) for i in range(len(fprs))], axis=0)
mean_tpr[-1] = 1.0
plt.plot(mean_fpr, mean_tpr, color='black', lw=2, linestyle='--', label=f'Mean ROC curve (AUC = {mean_auc:.2f})')

plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")

# Save the plot as a TIFF file with high quality
plt.savefig('roc_curve_cv_SVM.tif', dpi=300, bbox_inches='tight')

# Show plot
plt.show()
#------------------------------------------------------------------------------




svm_model.fit(X_train_scaled_enet_selected, y_train) #(Chnage the model seprate figure)
    
    # Predict and evaluate
y_pred = svm_model.predict(X_test_scaled_enet_selected)
y_proba_ev = svm_model.predict_proba(X_test_scaled_enet_selected)[:, 1]
    
accuracy = accuracy_score(y_test, y_pred)
auc = roc_auc_score(y_test, y_proba_ev)
results[name] = (accuracy, auc)
#print(f"{name} - Accuracy: {acc:.2f}, AUC: {auc:.2f}")

fpr, tpr, _ = roc_curve(y_test, y_proba_ev)
plt.figure()
plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (area = {auc:.2f})')
plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")



# Save the plot with high quality
plt.savefig('roc_curve_SVM.tif', dpi=300, bbox_inches='tight')
# Show plot
plt.show()



#------------------------------------------------------------------------------
# Save the trained SVM model
joblib.dump(svm_model, 'svm_model.pkl')

# Initialize and fit the scaler
scaler = StandardScaler().fit(X_train_scaled_enet_selected)

# Save the scaler
joblib.dump(scaler, 'scaler.pkl')


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import (
    accuracy_score, mean_squared_error, f1_score, precision_score,
    recall_score, confusion_matrix, matthews_corrcoef, roc_auc_score, roc_curve
)

# ---- Bootstrap CI function ----
def bootstrap_ci(y_true, y_pred, y_proba, n_bootstrap=1000, alpha=0.95):
    metrics = {
        "Accuracy": lambda yt, yp, yp_prob: accuracy_score(yt, yp),
        "AUC": lambda yt, yp, yp_prob: roc_auc_score(yt, yp_prob),
        "F1": lambda yt, yp, yp_prob: f1_score(yt, yp),
        "Precision": lambda yt, yp, yp_prob: precision_score(yt, yp),
        "Recall": lambda yt, yp, yp_prob: recall_score(yt, yp),
        "MCC": lambda yt, yp, yp_prob: matthews_corrcoef(yt, yp),
    }
    
    ci = {}
    n = len(y_true)
    for name, func in metrics.items():
        bootstrapped_scores = []
        for _ in range(n_bootstrap):
            indices = np.random.choice(range(n), n, replace=True)
            score = func(y_true[indices], y_pred[indices], y_proba[indices])
            bootstrapped_scores.append(score)
        lower = np.percentile(bootstrapped_scores, (1-alpha)/2 * 100)
        upper = np.percentile(bootstrapped_scores, (1 + alpha)/2 * 100)
        ci[name] = (lower, upper)
    return ci

# ---- Train & evaluate models ----
results = {}
for name, model in models.items():
    model.fit(X_train_scaled_enet_selected, y_train)

    # Predictions
    y_pred = model.predict(X_test_scaled_enet_selected)
    y_proba = model.predict_proba(X_test_scaled_enet_selected)[:, 1]

    # Metrics
    accuracy = accuracy_score(y_test, y_pred)
    mse = mean_squared_error(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    cm = confusion_matrix(y_test, y_pred)
    tn, fp, fn, tp = cm.ravel()
    sensitivity = recall
    specificity = tn / (tn + fp)
    mcc = matthews_corrcoef(y_test, y_pred)
    auc = roc_auc_score(y_test, y_proba)

    # Bootstrap CI
    ci_dict = bootstrap_ci(np.array(y_test), np.array(y_pred), np.array(y_proba), n_bootstrap=1000, alpha=0.95)

    results[name] = {
        "Accuracy": f"{accuracy:.3f} [{ci_dict['Accuracy'][0]:.3f}-{ci_dict['Accuracy'][1]:.3f}]",
        "AUC": f"{auc:.3f} [{ci_dict['AUC'][0]:.3f}-{ci_dict['AUC'][1]:.3f}]",
        "F1 Score": f"{f1:.3f} [{ci_dict['F1'][0]:.3f}-{ci_dict['F1'][1]:.3f}]",
        "Precision": f"{precision:.3f} [{ci_dict['Precision'][0]:.3f}-{ci_dict['Precision'][1]:.3f}]",
        "Recall": f"{recall:.3f} [{ci_dict['Recall'][0]:.3f}-{ci_dict['Recall'][1]:.3f}]",
        "Sensitivity": sensitivity,
        "Specificity": specificity,
        "MCC": f"{mcc:.3f} [{ci_dict['MCC'][0]:.3f}-{ci_dict['MCC'][1]:.3f}]",
        "MSE": mse
    }

# ---- Print results ----
for name, metrics in results.items():
    print(f"{name}:")
    for metric, value in metrics.items():
        print(f"  {metric}: {value}")

# # ---- ROC with bootstrap CI ----
n_bootstrap = 1000
rng = np.random.RandomState(42)

for name, model in models.items():
    # Get prediction probabilities
    y_proba = model.predict_proba(X_test_scaled_enet_selected)[:, 1]
    auc = roc_auc_score(y_test, y_proba)

    # Compute base ROC
    fpr, tpr, _ = roc_curve(y_test, y_proba)

    # Bootstrap TPRs
    bootstrapped_tprs = []
    for i in range(n_bootstrap):
        indices = rng.randint(0, len(y_proba), len(y_proba))
        if len(np.unique(y_test[indices])) < 2:  # skip if only one class
            continue
        fpr_i, tpr_i, _ = roc_curve(y_test[indices], y_proba[indices])
        bootstrapped_tprs.append(np.interp(fpr, fpr_i, tpr_i))

    tprs = np.array(bootstrapped_tprs)
    mean_tpr = tprs.mean(axis=0)
    std_tpr = tprs.std(axis=0)
    tpr_upper = np.minimum(mean_tpr + 1.96 * std_tpr, 1)
    tpr_lower = np.maximum(mean_tpr - 1.96 * std_tpr, 0)

    # ---- Plot for each model ----
    plt.figure()
    plt.plot(fpr, tpr, color="blue", lw=2, label=f"ROC curve (AUC = {auc:.2f})")
    plt.fill_between(fpr, tpr_lower, tpr_upper, color="blue", alpha=0.2, label="95% CI")
    plt.plot([0, 1], [0, 1], color="gray", lw=2, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve with 95% CI - {name}")
    plt.legend(loc="lower right")

    # save each model’s ROC+CI separately
    plt.savefig(f"roc_curve_with_CI_{name}.tif", dpi=1000, bbox_inches="tight")
    plt.close()











