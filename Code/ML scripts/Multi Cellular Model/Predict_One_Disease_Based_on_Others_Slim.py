import pandas as pd
import os
from sklearn.ensemble import RandomForestClassifier

"---------------------------- Load Data --------------------------------------------"

path= os.path.join('..', '..', 'Output', 'Full_scRNA_Datasets', 'Transfer_scRNA_Dataset_Fraction_Residuals.csv')
Full_scRNA = pd.read_csv(path)
print(Full_scRNA)

path= os.path.join('..', '..', 'Output', 'Runing_ML_scRNA_Datasets', 'Residuals_0.975_Transfer_scRNA_Synthetic_Cell_Proportion_10_Fraction.csv')
Synthetic_scRNA = pd.read_csv(path)
print(Synthetic_scRNA)

cols = list(Synthetic_scRNA)
print(cols)
relevant_features_dict ={'expression':'Expression', 'Fraction':'Fraction', 'preferential_expression':'Preferential expression', 'tipa_max':'Process activity score (relative)', 'num_interactors_dif_med':'Num. PPIs (relative)', 'diff_net_max': 'Differential PPI', 'num_interactors': 'Num. PPIs',  'paralogs_ratio_highest_identity':'Paralogs ratio'}

non_relevant_cols = ['Gene_ID', 'Cell_Type', 'Max_Standardized_Residuals', 'Affected_Cell']
relevant_cols = [c for c in cols if (c not in non_relevant_cols) & (c in relevant_features_dict)]
print('relevant_cols: ', len(relevant_cols), relevant_cols)

path= os.path.join('..', '..', 'Data', 'Diseases', 'Edited_OMIM_Morbidmap_Phenotypic_Series.csv')
PS_OMIM_Data = pd.read_csv(path, index_col=0)
print(PS_OMIM_Data)

" ----------------------------- Data Engineering ---------------------------"

Full_scRNA.fillna(0, inplace=True)
Synthetic_scRNA.fillna(0, inplace=True)

" ------------------ Predict Disease Genes in Cellular context ------------- "

disease_phenotypes = PS_OMIM_Data['# Phenotypic Series Number'].unique()

for ps in disease_phenotypes:

    ps_name = PS_OMIM_Data['Main_Phenotype'][PS_OMIM_Data['# Phenotypic Series Number'] == ps].tolist()[0]

    ps_name = ps_name.replace('?', '')
    ps_name = ps_name.replace('-', ' ')
    ps_name = ps_name.replace(',', ' ')
    ps_name = ps_name.replace("/", ' ')
    ps_name = ps_name.replace('  ', ' ')

    print(ps_name, ps)
    ps_genes = PS_OMIM_Data['Ensembl Gene ID (Ensembl)'][(PS_OMIM_Data['# Phenotypic Series Number'] == ps) & (PS_OMIM_Data['Certainty'] == 'High')].tolist()
    print('ps_genes', len(ps_genes), ps_genes)

    if len(ps_genes) > 0:

        Train_set = Synthetic_scRNA[~Synthetic_scRNA['Gene_ID'].isin(ps_genes)].copy(deep=True)
        Test_set = Full_scRNA[Full_scRNA['Gene_ID'].isin(ps_genes)].copy(deep=True)

        X_train = Train_set[relevant_cols]
        y_train = Train_set['Affected_Cell']
        X_test = Test_set[relevant_cols]

        model = RandomForestClassifier(random_state=1234)  # select the model
        model.fit(X_train, y_train)  # train the model
        y_pred = model.predict(X_test)  # predict the test data
        predictions_proba = model.predict_proba(X_test)
        pred_true = predictions_proba[:, 1]
        Test_set['scRNA_ML_Probability'] = pred_true

        try:
            path = os.path.join('..', '..', 'Output', 'Transfer_Learning_Diseases', 'Prediction_Slim', ps_name + '_' + ps + '_' + str(len(ps_genes)) + '_Cell_Specific_Probabilities_Slim.csv')
            Test_set.to_csv(path, index=None)

        except:
            print('$$$ Problem in : ', ps_name)