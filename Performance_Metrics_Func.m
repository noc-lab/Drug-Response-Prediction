%******************************************************************************************
%******************************************************************************************
% Matlab implementation of Drug response prediction algorithm introduced in the paper 
% "Machine Learning for Pharmacogenomics and Personalized Medicine: A Ranking Model for
% Drug Sensitivity Prediction" by Shahabeddin Sotudian  and Ioannis Ch.Paschalidis.
%******************************************************************************************

function [CI,SCI,AHat5,AHat10] = Performance_Metrics_Func(True_Rank,Predicted_Rank,sesetive_drugs)
% concordance index
CI=concordance_index(True_Rank,Predicted_Rank);

% concordance index for sensetives
Index_of_sensetives=find(sesetive_drugs);
[~, ~, sss1] = unique(True_Rank(:,Index_of_sensetives));
[~, ~, sss2] = unique(Predicted_Rank(:,Index_of_sensetives));
SCI=concordance_index(sss1',sss2');

% Average Hit @ K    
AHat5=   Average_Hit_AT_K(True_Rank,Predicted_Rank,5,sesetive_drugs);
AHat10=   Average_Hit_AT_K(True_Rank,Predicted_Rank,10,sesetive_drugs);

% Functions
function [ci] = concordance_index(targets, predicts)
  [~,i] = sort(targets,'descend');
  predicts = predicts(i);
  n = length(targets);
  total  = 0;
  norm_z = 0;
  for j=1:n
      for k=j:n
          if j ~= k
              h = step_function(predicts(j) - predicts(k));
              total = total + h;
              norm_z = norm_z + 1;
          end
      end
  end
  ci = total / norm_z;
end

 function h = step_function(diff)
    if diff > 0
        h = 1;
    elseif diff == 0
        h = 0.5;
    else
        h = 0;
    end
 end

function [OUT] = Average_Hit_AT_K(True_Rank,Predicted_Rank,K,sesetive_drugs)
    for i=1:K
    Index(i)=find(Predicted_Rank==i);
    end

    OUT=sum(sesetive_drugs(:,Index));

end
end