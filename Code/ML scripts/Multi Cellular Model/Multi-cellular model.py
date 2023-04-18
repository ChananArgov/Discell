
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
import pickle



" -------------------------------------- Load Data ------------------------------"

path= os.path.join('..', '..', 'Data', 'Residuals_0.975_Transfer_scRNA_Synthetic_Cell_Proportion_10_Fraction.csv')
Synthetic_Dataset_Residual = pd.read_csv(path)
print(Synthetic_Dataset_Residual)

"--------------------------------- Rename Columns  ---------------------------"

rename_dict = {'expression': 'Expression', 'Fraction':'Fraction',
               'tipa_max':'Biological processes', 'paralogs_ratio_highest_identity':'Paralogs compensation', 
               'preferential_expression':'Preferential expression', 'diff_net_max': 'Differential interactors', 
               'num_interactors':'PPIs'}

Synthetic_Dataset_Residual.rename(columns=rename_dict, inplace=True)

"--------------------------------- Process Data ---------------------------"

relevant_cols = [c for c in rename_dict.values() if c not in ['Gene_ID', 'Affected_Cell', 'Max_Standardized_Residuals', 'Cell_Type']]
print('relevant_cols', len(relevant_cols), relevant_cols)

Synthetic_Dataset_Residual[relevant_cols] = Synthetic_Dataset_Residual[relevant_cols].apply(pd.to_numeric, errors='coerce')
Synthetic_Dataset_Residual[relevant_cols] = Synthetic_Dataset_Residual[relevant_cols].fillna(Synthetic_Dataset_Residual[relevant_cols].median())


"--------------------------------- Cell Groups ---------------------------"

cell_types = Synthetic_Dataset_Residual['Cell_Type'][Synthetic_Dataset_Residual['Affected_Cell'] == True].unique()
print('cell_types', len(cell_types), cell_types)
general_cells = set([c.split('-')[1].strip().replace('cells', '') for c in cell_types])
general_cells = [c for c in general_cells if 'neurons' not in c]
general_cells.append('neurons')
print('general_cells', len(general_cells), general_cells)


precentage = 9
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 12))
results_dict = {'Cell_Type':[], 'ROC_AUC':[], 'PR_AUC':[], 'Disease_Causing_Genes':[]}


for gc in general_cells:

    print('@', gc)
    Test_set = Synthetic_Dataset_Residual[Synthetic_Dataset_Residual['Cell_Type'].str.contains(gc)]
    Train_set = Synthetic_Dataset_Residual[~Synthetic_Dataset_Residual['Cell_Type'].str.contains(gc)]
    test_positive = len(Test_set[Test_set['Affected_Cell']==True])
    X_train = Train_set[relevant_cols]
    y_train = Train_set['Affected_Cell']
    X_test = Test_set[relevant_cols]
    y_test = Test_set['Affected_Cell']

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
    pr_auc1 = auc(recall, prec)

    ax1.plot(fpr, tpr, label='%s ROC_AUC (area = %0.2f), Genes:%s' % (gc, roc_auc1, test_positive))
    ax2.plot(recall, prec, label='%s PR_AUC = %0.2f, Genes:%s' % (gc, pr_auc1, test_positive))

    results_dict['Cell_Type'].append(gc)
    results_dict['ROC_AUC'].append(roc_auc1)
    results_dict['PR_AUC'].append(pr_auc1)
    results_dict['Disease_Causing_Genes'].append(test_positive)


print(results_dict)
Results = pd.DataFrame.from_dict(results_dict)
print(Results)
path = os.path.join('..', '..', 'Results', 'Multi_Cellular_Model_LOOCV_AUC.csv')
Results.to_csv(path, index=None)

ax1.plot([0, 1], [0, 1], 'r--')
ax1.set_xlabel('1-Specificity(False Positive Rate)')
ax1.set_ylabel('Sensitivity(True Positive Rate)')
ax1.set_title('Receiver Operating Characteristic')
ax1.legend(loc="lower right", fontsize='small')

ax2.set_xlabel('Recall')
ax2.set_ylabel('Precision')
ax2.set_title('Precision-Recall curve')
ax2.axhline(y=1 / (precentage + 1), color='red', linestyle='--',
            label=r'Causal genes frequency = %0.2f' % (1 / (precentage + 1)))

ax2.legend(loc="lower right", fontsize='small')
plt.suptitle('General scRNA Model, Leave One Cell Out Residuals')

path = os.path.join('..', '..', 'Results',
                    'Multi_Cellular_Model_LOOCV_AUC_Fig.pdf')
plt.savefig(path)
plt.close()


"----------------------------- Model Explanation With SHAP ---------------------------"

X = Synthetic_Dataset_Residual[relevant_cols]
y = Synthetic_Dataset_Residual['Affected_Cell']

model = RandomForestClassifier(random_state=1234)  # select the model
model.fit(X, y)  # train the model

explainer = shap.TreeExplainer(model)

shap_values = explainer.shap_values(X)

shap.summary_plot(shap_values[1], X, show=False, plot_type='bar')

plt.tight_layout()
path = os.path.join('..', '..', 'Results', 'Multi_Cellular_Model_SHAP_Summery_Plot_bar.pdf')
plt.savefig(path)
plt.close()



