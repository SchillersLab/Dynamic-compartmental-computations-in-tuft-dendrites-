f = readtable('\\jackie-analysis\e\Shay\4482\2019-02-28_4482_freerun_aligned\Analysis\N2\Structural_VS_Functional\final\Run1\no_behave\Pearson\SP\cluster4\ByH\BetweenAndWithinSubTrees\AVSDForROI_1_Depth2_bp230.csv')
g = readtable('\\jackie-analysis\e\Shay\4482\2019-02-28_4482_freerun_aligned\Analysis\N2\Structural_VS_Functional\final\Run1\no_behave\Pearson\SP\cluster4\ByH\BetweenAndWithinSubTrees\AVSDForROI_3_Depth2_bp863.csv');
test1 = [g.Var1; f.Var1];
test2 = [g.Var2; f.Var2];
fitglm(test1, test2);
t = ans.Rsquared.Ordinary
yfit = predict(ans, test1);
plot(test1, yfit)