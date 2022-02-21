function SVMForAllData()
    mantelRunnerLocation = '\\jackie-analysis\e\Shay\RunnersLocationSummary.xlsx';
    sheetName = 'RunOnlyTuftPrecentageType1';
    tableResults1 = readtable(mantelRunnerLocation,'Sheet',sheetName);
    sheetName = 'RunOnlyTuftPrecentageType2';
    tableResults2 = readtable(mantelRunnerLocation,'Sheet',sheetName);

    outputfolder = '\\jackie-analysis\e\Shay\StatisticSummary\ByH_SVM\';

    for i = 1:size(tableResults1, 1)  
        for ci = 0:4
            load([tableResults1.RunLocation{i}, '\roiActivityRawData.mat'], 'selectedROISplitDepth1');
            load([tableResults1.RunLocation{i}, sprintf('\\roiActivityRaw_ByH_ByEvents_%d.mat', ci)], 'tmp');

            tmp = tmp';
            [~, ACC2D_depth1] = evalc("svmClassifyAndRand(tmp, selectedROISplitDepth1, selectedROISplitDepth1, 10, '', 1, 0)");
            chanceCalc = hist(selectedROISplitDepth1, unique(selectedROISplitDepth1));
            chanceCalc = chanceCalc/sum(chanceCalc);
            tableResults1.chanceLevelSVM(i) = max(chanceCalc);
            tableResults1.(sprintf('SVMAccuracy_%d', ci))(i) = ACC2D_depth1.mean;
            tmp = [];
            selectedROISplitDepth1 = [];
        end
    end

    for i = 1:size(tableResults2, 1)  
        for ci = 0:4
            load([tableResults2.RunLocation{i}, '\roiActivityRawData.mat'], 'selectedROISplitDepth1');
            load([tableResults2.RunLocation{i}, sprintf('\\roiActivityRaw_ByH_ByEvents_%d.mat', ci)], 'tmp');

            tmp = tmp';
            [~, ACC2D_depth1] = evalc("svmClassifyAndRand(tmp, selectedROISplitDepth1, selectedROISplitDepth1, 10, '', 1, 0)");
            chanceCalc = hist(selectedROISplitDepth1, unique(selectedROISplitDepth1));
            chanceCalc = chanceCalc/sum(chanceCalc);
            tableResults2.chanceLevelSVM(i) = max(chanceCalc);
            tableResults2.(sprintf('SVMAccuracy_%d', ci))(i) = ACC2D_depth1.mean;
            tmp = [];
            selectedROISplitDepth1 = [];
        end
    end
end