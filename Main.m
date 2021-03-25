%******************************************************************************************
%******************************************************************************************
% Matlab implementation of Drug response prediction algorithm introduced in the paper 
% "Machine Learning for Pharmacogenomics and Personalized Medicine: A Ranking Model for
% Drug Sensitivity Prediction" by Shahabeddin Sotudian  and Ioannis Ch.Paschalidis.
% -----------------------------------------------------------------------------------------
% Usage

% Inputs:
%   CellDrug_Data:      Drug sensitivity matrix where Rows are cell lines, Columns are drugs
%
%   CellGene_Data:      Gene expression matrix where each row represents gene expression of a cell line
%
%   Threshold:          Threshold to define sensitive drugs.
%
%   NonZeroElements:    After the data augmentation step, we can use all the elements or we can select
%                       a subset of indices to reduce the complexity of the model and improve its performance. 
%                       This can be tuned using cross-validation or recursive feature elimination approaches.
%
%   Lambda              Please refer to Equation 6 in the manuscript.
%
%   Phi_1               Please refer to Equation 1 in the manuscript.
%
%   Phi_2               Please refer to Equation 1 in the manuscript.

%******************************************************************************************

clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data and Parameters --------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
CellDrug_Data= % Load CellDrug Data
CellGene_Data= % Load CellGene Data
Number_Drugs=size(CellDrug_Data,2);   
Number_Genes=size(CellGene_Data,2);  

%% Parameters ---------------------------------------------------------
Threshold=20;
NonZeroElements = (1:(Number_Genes+Number_Drugs)^2)';  
Lambda=[1,10];
Phi_1=[0.1, 0.3];
Phi_2=[0.1, 0.3];

% Train Data
    CellDrug=100*CellDrug_Data(1:round(size(CellDrug_Data,1)*0.8),1:Number_Drugs);                
    CellGene=CellGene_Data(1:round(size(CellDrug_Data,1)*0.8),1:Number_Genes);

% Validation DATA 
    CellDrug_Valid=100*CellDrug_Data((round(size(CellDrug_Data,1)*0.8)+1):end,1:Number_Drugs);
    CellGene_Valid=CellGene_Data((round(size(CellDrug_Data,1)*0.8)+1):end,1:Number_Genes);
                  
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Training --------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Table_Final_Results] = Drug_Response_Prediction_Func(CellDrug,CellGene,CellDrug_Valid,CellGene_Valid,NonZeroElements,Phi_1,Phi_2,Threshold,Lambda);
disp('**************************************************************************************************')
disp('Drug Response Prediction - Model Performance:')
disp('**************************************************************************************************')
disp(Table_Final_Results)
disp('**************************************************************************************************')
           

    
    

