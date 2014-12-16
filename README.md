GIREN
=====

This directory provides the matlab code for GIREN approach, which integrate gene measurements and gene-gene interaction information under an elastic net formulation.

Please contact Dr. Lin Zhang <lin.zhang@cumt.edu.cn> if you have any questions regarding this package. 


For more technical details about the GIREN method, please refer to: 
 
Cancer Progression Prediction Using Gene Interaction Regularized Elastic Net 
Lin Zhang, Hui Liu, Yufei Huang, Xuesong Wang, Yidong Chen and Jia Meng*

Abstract
Different types of genomic aberration may simultaneously contribute to tumorigenesis. To obtain a more accurate prognostic assessment to guide therapeutic regimen choice for cancer patients, the heterogeneous multi-omic data should be integrated harmoniously, which can often be difficult. For this purpose, we propose a Gene Interaction Regularized Elastic Net (GIREN) model that predicts clinical outcome by integrating multiple data types. GIREN conveniently embraces both gene measurements and gene-gene interaction information under an elastic net formulation, enforcing structure sparsity and the “grouping effect” in solution to select the discriminant features with prognostic value. An iterative gradient descent algorithm is also developed to solve the model with regularized optimization. GIREN was applied to human ovarian cancer and breast cancer datasets obtained from The Cancer Genome Atlas, respectively. Result shows that, the proposed GIREN algorithm consistently outperforms competing algorithms (LASSO, Elastic Net and superPC) in predicting cancer progression on both two datasets under different settings, suggesting a promising direction for more effective integration of gene measurement and interaction information.
