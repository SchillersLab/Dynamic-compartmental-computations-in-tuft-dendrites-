
valsAll=zeros(5,40,52);
matR2 = [];matR1=[];
for i = 1:size(tableResults1, 1)
for ci = 0:4
load([tableResults1.RunLocation{i}, '\roiActivityRawData.mat'], 'selectedROISplitDepth1');
load([tableResults1.RunLocation{i}, sprintf('\\roiActivityRaw_ByH_ByEvents_%d.mat', ci)], 'tmp');
[pcares.embedding, ~, vals] = pca(tmp);
pcares.eigs = vals;

valsDiff = diff(vals./sum(vals));
valsCumsum = cumsum(vals./sum(vals));

% indexT = find(abs(valsDiff) <= 0.05,1);
indexTsec = find(valsCumsum >= 0.90, 1);

valsAll(ci+1, 1:length(vals), i) = vals./sum(vals);   

% if indexT == length(vals)
%     error('bag');
% end

% if indexT < indexTsec
    indexT = indexTsec;
% end

if indexT == 1
    indexT = 2;
end

matR1(i, ci+1) = indexT; 

embedding = pcares.embedding(:,1:indexT);
% embedding = pcares.embedding(:,1:2);
[~, ACC2D_depth1] = evalc("svmClassifyAndRand(embedding, selectedROISplitDepth1, selectedROISplitDepth1, 10, '', 1, 0)");
chanceCalc = hist(selectedROISplitDepth1, unique(selectedROISplitDepth1));
chanceCalc = chanceCalc/sum(chanceCalc);
tableResults1.chanceLevelSVM(i) = max(chanceCalc);
tableResults1.(sprintf('SVMAccuracy_%d', ci))(i) = ACC2D_depth1.mean;
tmp = [];tmp2 = [];
selectedROISplitDepth1 = [];
end
end

for i = 1:size(tableResults2, 1)
for ci = 0:4
load([tableResults2.RunLocation{i}, '\roiActivityRawData.mat'], 'selectedROISplitDepth1');
load([tableResults2.RunLocation{i}, sprintf('\\roiActivityRaw_ByH_ByEvents_%d.mat', ci)], 'tmp');
[pcares.embedding, ~, vals] = pca(tmp);
pcares.eigs = vals;

valsDiff = diff(vals./sum(vals));
valsCumsum = cumsum(vals./sum(vals));

% indexT = find(abs(valsDiff) <= 0.05,1);
indexTsec = find(valsCumsum >= 0.90, 1);

valsAll(ci+1, 1:length(vals), i+27) = vals./sum(vals);

% if indexT == length(vals)
%     error('bag');
% end


% if indexT < indexTsec
    indexT = indexTsec;
% end

if indexT == 1
    indexT = 2;
end
matR2(i, ci+1) = indexT; 

embedding = pcares.embedding(:,1:indexT);
% embedding = pcares.embedding(:,1:2);
[~, ACC2D_depth1] = evalc("svmClassifyAndRand(embedding, selectedROISplitDepth1, selectedROISplitDepth1, 10, '', 1, 0)");
chanceCalc = hist(selectedROISplitDepth1, unique(selectedROISplitDepth1));
chanceCalc = chanceCalc/sum(chanceCalc);
tableResults2.chanceLevelSVM(i) = max(chanceCalc);
tableResults2.(sprintf('SVMAccuracy_%d', ci))(i) = ACC2D_depth1.mean;
tmp = [];tmp2 = [];
selectedROISplitDepth1 = [];
end
end

valsAll2 = cumsum(valsAll, 2);
for ci = 0:4
    figure; hold on;
    errorbar(mean(valsAll(ci+1, :, :), 3), std(valsAll(ci+1, :, :),0, 3), 'DisplayName', sprintf('cluster %d', ci))
    ylim([0,1])
    xlim([0,10])
    xlim([0,20])
    figure; hold on;
    errorbar(mean(valsAll2(ci+1, :, :), 3), std(valsAll2(ci+1, :, :),0, 3), 'DisplayName', sprintf('cluster %d', ci))
    ylim([0,1])
    xlim([0,10])
    xlim([0,20])
end

pcaMeanAll(1) = {[tableResults1.SVMAccuracy_0, tableResults1.SVMAccuracy_1,tableResults1.SVMAccuracy_2,tableResults1.SVMAccuracy_3,tableResults1.SVMAccuracy_4]};
pcaMeanAll(2) = {[tableResults2.SVMAccuracy_0, tableResults2.SVMAccuracy_1,tableResults2.SVMAccuracy_2,tableResults2.SVMAccuracy_3,tableResults2.SVMAccuracy_4]};
pcaChance(1) = {mean(tableResults1.chanceLevelSVM)};
pcaChance(2) = {mean(tableResults2.chanceLevelSVM)};
f = figure;hold on;
title('PCA All accuracy Mean')
groupsBox1 = {};
groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-All'};
groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster1'};
groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster2'};
groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster3'};
groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster4'};
groupsBox2 = {};
groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-All'};
groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster1'};
groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster2'};
groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster3'};
groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster4'};
meanAll(1:5) = mean(pcaMeanAll{1});
meanAll(6:10) = mean(pcaMeanAll{2});
chaAll(1:5) = ones(1,5)*(pcaChance{1});
chaAll(6:10) = ones(1,5)*(pcaChance{2});
bC = boxchart(categorical(groupsBox1),[pcaMeanAll{1}(:,1); pcaMeanAll{1}(:,2); pcaMeanAll{1}(:,3); pcaMeanAll{1}(:,4); pcaMeanAll{1}(:,5)]);
bC2 = boxchart(categorical(groupsBox2),[pcaMeanAll{2}(:,1); pcaMeanAll{2}(:,2); pcaMeanAll{2}(:,3); pcaMeanAll{2}(:,4); pcaMeanAll{2}(:,5)]);
bC.BoxFaceColor = [0,0,0];
bC.BoxFaceAlpha = 0.4;
bC.MarkerColor = [0,0,0];
bC2.BoxFaceColor = [0,0,255] ./ 255;
bC2.MarkerColor = [0,0,255] ./ 255;
bC2.BoxFaceAlpha = 0.4;
plot(meanAll, '-*k');
plot(chaAll, '--k')
ylim([0,1]);

[p1, h1] = ranksum(pcaMeanAll{1}(:,1),  pcaMeanAll{2}(:,1));
[p2, h2] = ranksum(pcaMeanAll{1}(:,2),  pcaMeanAll{2}(:,2));
[p3, h3] = ranksum(pcaMeanAll{1}(:,3),  pcaMeanAll{2}(:,3));
[p4, h4] = ranksum(pcaMeanAll{1}(:,4), pcaMeanAll{2}(:,4));
[p5, h5] = ranksum(pcaMeanAll{1}(:,5), pcaMeanAll{2}(:,5));

text = '';
text = strcat(text, sprintf('Cluster all h1 %f, p %f \\n', h1,p1));
text = strcat(text, sprintf('Cluster2 h1 %f, p %f \\n', h2,p2));
text = strcat(text, sprintf('Cluster3 h1 %f, p %f \\n', h3,p3));
text = strcat(text, sprintf('Cluster4 h1 %f, p %f \\n', h4,p4));
text = strcat(text, sprintf('Cluster4 h1 %f, p %f \\n', h5,p5));

fid=fopen(fullfile('\\jackie-analysis\E\Shay\StatisticSummary\ByH\SVM', 'statisticranksumPCASmallVsBig.txt'),'w');
fprintf(fid, text);
fclose(fid);

textAnova = '';
[p,~,statsM] = anova1(pcaMeanAll{1}, {'All', 'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4'}, 'off');
if p < 0.05
    [c,~,~,groupnames] = multcompare(statsM, 'Display', 'off');
    for j = 1:size(c, 1)
        textAnova = strcat(textAnova, sprintf('%s vs %s - pValue : %.4f, mean Diff : %.4f, CI_L : %.4f, CI_H : %.4f \\n ', ...
            groupnames{c(j, 1)}, groupnames{c(j, 2)}, c(j, 6), c(j, 4), c(j, 3), c(j, 5)));
    end
else
    textAnova = 'Anova not pass';
end

fid=fopen(fullfile('\\jackie-analysis\E\Shay\StatisticSummary\ByH\SVM', ['statistic_Big.txt']),'w');
fprintf(fid, textAnova);
fclose(fid);

textAnova = '';
[p,~,statsM] = anova1(pcaMeanAll{2}, {'All', 'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4'}, 'off');
if p < 0.05
    [c,~,~,groupnames] = multcompare(statsM, 'Display', 'off');
    for j = 1:size(c, 1)
        textAnova = strcat(textAnova, sprintf('%s vs %s - pValue : %.4f, mean Diff : %.4f, CI_L : %.4f, CI_H : %.4f \\n ', ...
            groupnames{c(j, 1)}, groupnames{c(j, 2)}, c(j, 6), c(j, 4), c(j, 3), c(j, 5)));
    end
else
    textAnova = 'Anova not pass';
end

fid=fopen(fullfile('\\jackie-analysis\E\Shay\StatisticSummary\ByH\SVM', ['statistic_Small.txt']),'w');
fprintf(fid, textAnova);
fclose(fid);