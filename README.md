# Discell
Discell is a ML-based approach aim to illuminate the cellular basis of hereditary diseases.

<img src="Concept Figure.png" alt="Concept Figure">

# Dataset
The dataset contains Discell features and gene labels per cell type, can be found [here](https://zenodo.org/record/5769155#.Yh9sEOhBwuU](https://zenodo.org/deposit/5769155).

<h2>Code</h2>
<table>
  <tr>
    <th>Script name</th>
    <th>Folder</th>
    <th>Contant</th>
    <th>Manuscript figures</th>
    <th>Comments</th>
  </tr>
  <tr>
    <td>Query_Specific_ML_Training_Server.py</td>
    <td>ML scripts/One Disease Model/Query_Specific_ML_Training_Server.py</td>
    <td>Disease specific ML model train, evaluation and explanation</td>
    <td>Fig. 2B-E</td>
    <td></td>
  </tr>  
  <tr>
    <td>Create_Query_Dataset_Server.py</td>
    <td>ML scripts\One Disease Model\Create_Query_Dataset_Server.py</td>
    <td>Creates the training dataset for the disease specific ML model</td>
    <td></td>
    <td>The input is list of relevant genes</td>
  </tr>
  <tr>
    <td>Investigate model predictions Autism.py</td>
    <td>Analysis Script</td>
    <td>ASD model further validation</td>
    <td>Fig. 2F</td>
    <td></td>
  </tr>
  <tr>
    <td>Slim_model_Performance_Plot.ipynb</td>
    <td>Analysis Script</td>
    <td>Multi-cellular model ROC-AUC, PR-AUC plot</td>
    <td>Fig. 3B</td>
    <td>Based on Table S4</td>
  </tr> 
  <tr>
    <td>Multi-cellular model.py</td>
    <td>ML scripts/Multi Cellular Model</td>
    <td>Multi-cellular model train, evaluation and explanation</td>
    <td>Fig. 3C</td>
    <td></td>
  </tr>
  <tr>
    <td>ML_Disease_Predictor_OMIM Slim Z score_Matrix.py</td>
    <td>Comparision</td>
    <td>Compare the multi-cellular model performance to the preferential based approach</td>
    <td>Fig. 3D</td>
    <td></td>
  </tr>
  <tr>
    <td>Diseases_Predictions_Heatmap_Slim.py</td>
    <td>Analysis Script</td>
    <td>Disease probabilities heatmap plot</td>
    <td>Fig. 4</td>
    <td></td>
  </tr>
  <tr>
    <td>Predict_One_Disease_Based_on_Others_Slim.py</td>
    <td>ML scripts/Multi Cellular Model</td>
    <td>Leave one disease out cross validation</td>
    <td>Refers to Fig. 5</td>
    <td></td>
  </tr>
    <tr>
    <td>Predict_One_Disease_Probabilities_Plot_Slim.py</td>
    <td>Analysis Script</td>
    <td>Disease genes probability across cell type boxplot</td>
    <td>Fig. 5 left panle</td>
    <td></td>
  </tr> 
  <tr>
    <td>Disease_Probability_UMAP_General.py</td>
    <td>Analysis Script</td>
    <td>Disease genes probability across cell type on UMAP</td>
    <td>Fig. 5 midel panle</td>
    <td></td>
  </tr> 
  <tr>
    <td>Give_MultiOmics_Net_Server.py</td>
    <td>Analysis Script</td>
    <td>Disease genes directed interactome</td>
    <td>Fig. 5 right panle</td>
    <td></td>
  </tr> 
</table>
</body>
</html>


# Cite
Please cite 'Illuminating hereditary disease mechanisms in cell-specific contexts through multi-omics data and machine learning. Argov et al, submitted.'

# Contact
Esti Yeger-Lotem, estiyl@bgu.ac.il
Chanan Argov, chanana@post.bgu.ac.il
