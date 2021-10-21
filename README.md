# Drug Response Prediction
Matlab implementation of drug response prediction algorithm introduced in the paper "Machine Learning for Pharmacogenomics and
Personalized Medicine: A Ranking Model for  Drug Sensitivity Prediction" by Shahabeddin Sotudian  and Ioannis Ch.Paschalidis.
It is infeasible to test many different chemotherapy drugs on actual patients in large clinical trials, which motivates computational methods with the ability to learn and exploit associations between drug effectiveness and patient characteristics. This work proposes a machine learning approach to infer robust predictors of drug responses from patient genomic information. Rather than predicting the exact drug response on a given cell line, we introduce an elastic-net regression methodology to compare a drug-cell line pair against an alternative pair. Using predicted pairwise comparisons we rank the effectiveness of different drugs on the same cell line. A total of 173 cell lines and 100 drug responses were used in various settings for training and testing the proposed models. By comparing our approach against twelve baseline methods, we demonstrate that it outperforms the state-of-the-art methods in the literature. In contrast to most other methods, the algorithm is able to maintain its high performance even when we use a large number of drugs and few cell lines.

 
 
 
## Citation

If you use the code, please cite this paper:

```text
@article{sotudian2021machine,
  title={Machine Learning for Pharmacogenomics and Personalized Medicine: A Ranking Model for Drug Sensitivity Prediction},
  author={Sotudian, Shahabeddin and Paschalidis, Ioannis CH},
  journal={IEEE/ACM Transactions on Computational Biology and Bioinformatics},
  year={2021},
  publisher={IEEE}
}
```
