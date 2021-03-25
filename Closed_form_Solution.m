%******************************************************************************************
%******************************************************************************************
% Matlab implementation of Drug response prediction algorithm introduced in the paper 
% "Machine Learning for Pharmacogenomics and Personalized Medicine: A Ranking Model for
% Drug Sensitivity Prediction" by Shahabeddin Sotudian  and Ioannis Ch.Paschalidis.
%******************************************************************************************
function [W_Tilde]= Closed_form_Solution(F,P,lambda,H)  
X_Tilde=cat(2,ones(size(F,1),1),F);
C2=eye(size(X_Tilde,2));    
C2(1,1)=0;          
%% Solution
W_Tilde=pinv((((X_Tilde')*X_Tilde)+(1/H)*lambda*C2))*(X_Tilde')*P;
end    