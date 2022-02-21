function mainRunnerActivitySetUp
    
    neuronTreePathSWC = "E:\Dropbox (Technion Dropbox)\Yara\Analysis\Yara's Data For Shay\SM01\06.27.19_FreeRun\SM01_06.27.19_FreeRunswcFiles\neuron_1.swc";
    
    outputpath = "E:\Dropbox (Technion Dropbox)\Yara\Analysis\Yara's Data For Shay\SM01\06.27.19_FreeRun\TreeCorrResultsNew\neuron_1\SetUP\";
    outputpath = char(outputpath);
   

    activityByCSV = false;
    neuronActiityPathCSV = 'D:\Shay\work\N\SMO1_ 15_8_19\15_08_19_aligned_Ch1_dffByRoi.csv';
    neuronActivityPathTPA = 'E:\Dropbox (Technion Dropbox)\Yara\Analysis\SM01\06.27.19_FreeRun\';
    
    threshold_std = 3;
    PeakWidth = 5;
    clusterPrecentage = 5;
    histBinWidth = 0.01;
    
    samplingRate = 20;
    
    [~, ~, selectedROITable] = loadSwcFile(neuronTreePathSWC, outputpath);
    
    selectedROI = selectedROITable.Name;
   
    
    if (activityByCSV)
        %     load roi activity file
        [roiActivity, roiActivityNames] = loadActivityFile(neuronActiityPathCSV, selectedROI);
    else
        [roiActivity, roiActivityNames] = loadActivityFileFromTPA(neuronActivityPathTPA, selectedROI, outputpath);
    end
        
    
%     Calc Activity Events Window
    [windowFULL, windowToPeak, locationPeaks, valuePeaks, activityClusterValue] = calcActivityEventsWindowsAndPeaks(roiActivity, outputpath, threshold_std, PeakWidth, clusterPrecentage, histBinWidth);
end