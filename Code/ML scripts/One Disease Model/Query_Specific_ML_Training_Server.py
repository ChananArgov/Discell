import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
import sklearn as skl
from sklearn.metrics import average_precision_score
import shap
import sys
import pickle

" -------------------------------------- Load Data ------------------------------"
file_name = sys.argv[1]
name = file_name.replace('.csv', '')
path= os.path.join('..', 'Output', name  + '_Query_Synthetic_Dataset.csv')

Query_Synthetic_Dataset = pd.read_csv(path)
print(list(Query_Synthetic_Dataset))

"--------------------------------- Relevant Features ---------------------------"
path= os.path.join('..',  'Data', 'Relevant_Features_Names_scRNA.csv')
Relevant_Features = pd.read_csv(path)
rename_dict = dict(zip(Relevant_Features['Feature'], Relevant_Features['Feature Name New']))

relevant_features = Relevant_Features['Feature Name New'][Relevant_Features['Relevant Feature'] == True].tolist()

"--------------------------------- Process Data --------------------------------"

Query_Synthetic_Dataset.rename(columns=rename_dict, inplace=True)
relevant_cols = [c for c in list(Query_Synthetic_Dataset) if c not in ['Gene_ID', 'Query_Gene'] and c in relevant_features]
print('relevant_cols', len(relevant_cols))

Query_Synthetic_Dataset[relevant_cols] = Query_Synthetic_Dataset[relevant_cols].apply(pd.to_numeric, errors='coerce')
Query_Synthetic_Dataset[relevant_cols] = Query_Synthetic_Dataset[relevant_cols].fillna(Query_Synthetic_Dataset[relevant_cols].median())


"-------------------------------- Run ML ---------------------------"

kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=1234)

X = Query_Synthetic_Dataset[relevant_cols]
y = Query_Synthetic_Dataset['Query_Gene']
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))
fold = 0
results_dict = {'Fold': [], 'ROC_AUC':[], 'PR_AUC':[]}

for train, test in kf.split(X, y):
    fold += 1
    train1 = train.tolist()
    test1 = test.tolist()
    X_train, X_test = X.iloc[train1], X.iloc[test1]

    y_train, y_test = y.iloc[train1], y.iloc[test1]

    model = RandomForestClassifier(random_state=1234)  # select the model
    model.fit(X_train, y_train)  # train the model
    y_pred = model.predict(X_test)  # predict the test data
    predictions_proba = model.predict_proba(X_test)
    clr = classification_report(y_test, y_pred, output_dict=True)
    precision = clr['True']['precision']
    recall_1 = clr['True']['recall']
    f1_score = clr['True']['f1-score']
    print('precision:', precision, 'recall: ', recall_1, 'f1_score: ', f1_score)
    print(model.classes_)
    roc_auc1 = roc_auc_score(y_test, predictions_proba[:, 1])
    pred_true = predictions_proba[:, 1]
    print('AUC = %0.2f' % roc_auc1)
    print('\n')
    fpr, tpr, _ = roc_curve(y_test, predictions_proba[:, 1], pos_label=model.classes_[1])
    prec, recall, _ = precision_recall_curve(y_test, predictions_proba[:, 1], pos_label=model.classes_[1])
    pr_auc1 = average_precision_score(y_test, pred_true)
    ax1.plot(fpr, tpr, label='Fold %s ROC_AUC (area = %0.2f)' % (fold, roc_auc1))
    ax2.plot(recall, prec, label='Fold %s PR_AUC = %0.2f' % (fold, pr_auc1))

    results_dict['Fold'].append(fold)
    results_dict['ROC_AUC'].append(roc_auc1)
    results_dict['PR_AUC'].append(pr_auc1)


print(results_dict)
Results = pd.DataFrame.from_dict(results_dict)
print(Results)

path = os.path.join('..', 'Output', name + '_Query_10_Fold_Performance.csv')
Results.to_csv(path, index=None)

ax1.plot([0, 1], [0, 1], 'r--')
ax1.set_xlabel('1-Specificity(False Positive Rate)')
ax1.set_ylabel('Sensitivity(True Positive Rate)')
ax1.set_title('Receiver Operating Characteristic')
ax1.legend(loc="lower right", fontsize='small')

ax2.set_xlabel('Recall')
ax2.set_ylabel('Precision')
ax2.set_title('Precision-Recall curve')
ax2.axhline(y=0.1, color='red', linestyle='--',
            label=r'Causal genes frequency = %0.2f' % (0.1))

ax2.legend(loc="lower right", fontsize='small')
plt.suptitle('scRNA ML Model for Query Genes')

path = os.path.join('..', 'Output', name + '_Query_10_Fold_Performance_AUC.jpg')
plt.savefig(path)
plt.close()

"------------------------------- SHAP Part --------------------------"

model = RandomForestClassifier(random_state=1234)  # select the model
model.fit(X, y)  # train the model

# Create object that can calculate shap values
explainer = shap.TreeExplainer(model)

path =  os.path.join('..', 'Output', name + '_Explainer.pkl')
with open(path, 'wb') as handle:
    pickle.dump(explainer, handle, protocol=pickle.HIGHEST_PROTOCOL)

# calculate shap values. This is what we will plot.
# Calculate shap_values for all of val_X rather than a single row, to have more data for plot.
shap_values = explainer.shap_values(X, check_additivity=False) # Worning!!!

# Make plot. Index of [1] is explained in text below.
shap.summary_plot(shap_values[1], X, show=False)

plt.tight_layout()
path = os.path.join('..', 'Output', name + '_Query_Model_SHAP_Summery_Plot.jpg')
plt.savefig(path)
