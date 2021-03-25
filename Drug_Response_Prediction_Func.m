
%******************************************************************************************
%******************************************************************************************
% Matlab implementation of Drug response prediction algorithm introduced in the paper 
% "Machine Learning for Pharmacogenomics and Personalized Medicine: A Ranking Model for
% Drug Sensitivity Prediction" by Shahabeddin Sotudian  and Ioannis Ch.Paschalidis.
% -----------------------------------------------------------------------------------------
% Usage

% Inputs:
%   CellDrug:      Drug sensitivity matrix where Rows are cell lines, Columns are drugs
%
%   CellGene:      Gene expression matrix where each row represents gene expression of a cell line
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
%
% Output
%   Performance_Table
%******************************************************************************************



function [Performance_Table] = Drug_Response_Prediction_Func(CellDrug,CellGene,CellDrug_Valid,CellGene_Valid,NonZeroElements,Phi_1,Phi_2,Thereshold,Lambda)

% Create main data X Y
N_drug=size(CellDrug,2); % Number of drugs
N_Cell=size(CellDrug,1); % Number of Cells
N_Gene=size(CellGene,2); % Number of Genes

X_int=cat(2,eye(N_drug),zeros(N_drug,N_Gene));

for i=1:N_Cell
        X2=X_int;
        Y2=zeros(N_drug,1);
          for j=1:N_drug
             X2(j,(N_drug+1):(N_drug+N_Gene))=CellGene(i,:);
             Y2(j,1)=CellDrug(i,j);
          end
     if i>1
           X=cat(1,X,X2);
           Y=cat(1,Y,Y2);
     else
           X=X2;  
           Y=Y2;
     end
end

% Pairwise comparison - the data augmentation step
N=N_Cell * N_drug;
Total_data=(N_Cell-1)*(N) +   (N_drug-1)*N;
Counter=1;
LL=1;
DD=1;
XX=zeros( Total_data ,  size(NonZeroElements,1) );
YY=zeros( Total_data ,  1 );

for i=1:N
    for j=1:N
        if i~=j
            if (X(j,1:N_drug)==X(i,1:N_drug)) 
               SpareVar=(reshape(  (X(i,:)')*X(j,:) ,[],1)');
               XX(Counter,:)=SpareVar(:,NonZeroElements);
               YY(Counter,:)=(Y(i)-Y(j));    
               zarib_ONE_index(LL)=Counter;
               LL=LL+1;
               Counter=Counter+1;

            end

            if (X(j,(N_drug+1):(N_Gene+N_drug))==X(i,(N_drug+1):(N_Gene+N_drug)))
               SpareVar=(reshape(  (X(i,:)')*X(j,:) ,[],1)');
               XX(Counter,:)=SpareVar(:,NonZeroElements);
               YY(Counter,:)=Y(i)-Y(j);    
               zarib_TWO_index(DD)=Counter;
               DD=DD+1;
               Counter=Counter+1;
            end
        end
    end 
end

Size_W=size(X,2);
clearvars -except Lambda XX YY NonZeroElements Size_W CellDrug_Valid CellGene_Valid X N_drug N_Gene Thereshold N_Cell N Phi_1 Phi_2 zarib_ONE_index zarib_TWO_index H Num_Nonzero_elements;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Learning --------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


counter=1;
for mm=1:size(Lambda,2)
    for zz=1:size(Phi_1,2)
        for rr=1:size(Phi_2,2)

            % Define Saay
            SAAY=(ones(size(YY)));
            SAAY(zarib_ONE_index)=Phi_1(zz);
            SAAY(zarib_TWO_index)=Phi_2(rr);

            % Define H
            H= (N_drug+N_Cell-2)*(N);  %  can use H=1 instead

            % Solve using closed form
            [W_Tilde]= Closed_form_Solution((XX),(YY.*SAAY),Lambda(mm),H);

            W=zeros(Size_W);
            W(NonZeroElements)=W_Tilde(2:size(W_Tilde,1));
            b=W_Tilde(1);

            % Validation  
            N_drug_Validation=size(CellDrug_Valid,2); % Number of drugs
            N_Cell_Validation=size(CellDrug_Valid,1); % Number of Cells
            N_Gene_Validation=size(CellGene_Valid,2); % Number of Genes

            for i=1:size(CellGene_Valid,1)
                    X_Validation=cat(2,eye(N_drug_Validation),repmat(CellGene_Valid(i,:),N_drug_Validation,1));
                    Score_Validation=zeros(1,size(X_Validation,1));
                    Y_Validation=Score_Validation;
                    Cat_X_XValidation=cat(1,X_Validation,X);
                    % SCORE
                            for i1=1:size(X_Validation,1)
                                        for j1=1:size(Cat_X_XValidation,1)
                                                if i1~=j1
                                                    if (X_Validation(i1,1:N_drug)==Cat_X_XValidation(j1,1:N_drug))      
                                                        Score_Validation(i1)=Score_Validation(i1)+(X_Validation(i1,:)*W)*(Cat_X_XValidation(j1,:)')+b; 
                                                    end
                                                    if ((X_Validation(i1,(N_drug+1):(N_Gene+N_drug))==Cat_X_XValidation(j1,(N_drug+1):(N_Gene+N_drug))))
                                                        Score_Validation(i1)=Score_Validation(i1)+(X_Validation(i1,:)*W)*(Cat_X_XValidation(j1,:)')+b; 
                                                    end

                                                end
                                        end

                                        Y_Validation(i1)=CellDrug_Valid(i,i1);
                            end

                    %  Score_Validation
                    Score_Validation_RAW(i,:)=Score_Validation;
                    [~,p12222] = sort(Score_Validation,'ascend');
                    r23333 = 1:length(Score_Validation);
                    ranking1_Validation=[];
                    ranking1_Validation(p12222) = r23333;
                    % Y_Validation
                    [~,p1] = sort(Y_Validation,'ascend');
                    r2 = 1:length(Y_Validation);
                    ranking2_Validation(p1) = r2;
                    Result_Validation=cat(1,ranking2_Validation,ranking1_Validation);
                    % Finding sensetive drugs
                    sesetive_drugs=CellDrug_Valid(i,:);
                    prctile(CellDrug_Valid(i,:),Thereshold);
                    ss1= find(CellDrug_Valid(i,:)< prctile(CellDrug_Valid(i,:),Thereshold));
                    sesetive_drugs(:,ss1)=1;
                    ss2= find(sesetive_drugs~=1);
                    sesetive_drugs(ss2)=0;
                    [CI(i),SCI(i),AHat5(i),AHat10(i)] = Performance_Metrics_Func(Result_Validation(1,:),Result_Validation(2,:),sesetive_drugs);
            end
            
            Selection_Results(counter,1)=counter;
            Selection_Results(counter,2)=Phi_1(zz);
            Selection_Results(counter,3)=Phi_2(rr);
            Selection_Results(counter,4)=Lambda(mm);
            Selection_Results(counter,5)=mean(CI);
            Selection_Results(counter,6)=mean(SCI);
            Selection_Results(counter,7)=mean(AHat5);
            Selection_Results(counter,8)=mean(AHat10);
            
            
            counter=counter+1;

        end
    end
end

%% Display results
Num = Selection_Results(:,1);
Phi_1 =Selection_Results(:,2);
Phi_2=Selection_Results(:,3);
Lambda = Selection_Results(:,4);
T_CI = Selection_Results(:,5);
T_SCI = Selection_Results(:,6);
T_AHat5 = Selection_Results(:,7);
T_AHat10 = Selection_Results(:,8);


Performance_Table = table(Num,Phi_1,Phi_2,Lambda,T_CI,T_SCI,T_AHat5,T_AHat10);

end


