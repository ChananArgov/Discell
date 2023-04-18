import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

" ------------------------- Load Data -------------------------- "

files_list = os.listdir(os.path.join('..', '..', 'Output', 'Transfer_Learning_Diseases', 'Prediction_Slim'))
print('files_list', len(files_list), files_list)

for file_name in files_list:

    short_name = file_name.replace('_Specific_Probabilities.csv', '')
    disease_name = file_name.split('_')[0]
    ps_num = file_name.split('_')[1]
    genes = file_name.split('_')[2]

    path = os.path.join('..', '..', 'Output', 'Transfer_Learning_Diseases', 'Prediction_Slim', file_name)
    Disease_Cell_Probabilities = pd.read_csv(path)

    print(Disease_Cell_Probabilities)
    plt.figure(0, figsize=(60, 12))


    Disease_Prediction_Mean = Disease_Cell_Probabilities.groupby('Cell_Type')['scRNA_ML_Probability'].mean()
    Disease_Prediction_Mean = pd.DataFrame(Disease_Prediction_Mean)
    Disease_Prediction_Mean['scRNA_ML_Probability'] = Disease_Prediction_Mean['scRNA_ML_Probability'] / Disease_Prediction_Mean['scRNA_ML_Probability'].max()
    Disease_Prediction_Mean.reset_index(inplace=True)
    Disease_Prediction_Mean['Color_Degree'] = 'silver'
    Disease_Prediction_Mean.loc[Disease_Prediction_Mean['scRNA_ML_Probability'] < 0.2, 'Color_Degree'] = 'lightgrey'
    Disease_Prediction_Mean.loc[Disease_Prediction_Mean['scRNA_ML_Probability'] > 0.6, 'Color_Degree'] = 'lightsalmon'
    Disease_Prediction_Mean.loc[Disease_Prediction_Mean['scRNA_ML_Probability'] > 0.9, 'Color_Degree'] = 'red'

    color_dict = {cell: Disease_Prediction_Mean[Disease_Prediction_Mean['Cell_Type'] == cell]['Color_Degree'].values[0] for cell in Disease_Cell_Probabilities['Cell_Type'].unique()}

    sns.boxplot(data=Disease_Cell_Probabilities, x='Cell_Type', y='scRNA_ML_Probability', palette=color_dict)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.title(disease_name + ' ' + ps_num + ' genes: ' + genes)
    try:
        path = os.path.join('..', '..', 'Output', 'Transfer_Learning_Diseases', 'Cell_Disease_Plots_Slim', short_name + '_Plot.png')
        plt.savefig(path)
        plt.close()
    except:
        print('*** Problem:', short_name)

    # break
