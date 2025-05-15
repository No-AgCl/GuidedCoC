function [TAB_Y, Eval_tab] = clu_eval(clu_Y_truth, clu_Y_bes)
sample = zeros(4,1);
sample(1)=Cal_Purity(clu_Y_truth, clu_Y_bes);
[sample(2),sample(3),~,~]=RandIndex(clu_Y_truth,clu_Y_bes);
sample(4) = Cal_NMI(clu_Y_truth, clu_Y_bes);
rowNames = {'Purity','ARI','RI','NMI'};
colNames = {'Y'};
Eval_tab = array2table(sample,'RowNames',rowNames,'VariableNames',colNames);

TAB_Y = crosstab(clu_Y_truth,clu_Y_bes);