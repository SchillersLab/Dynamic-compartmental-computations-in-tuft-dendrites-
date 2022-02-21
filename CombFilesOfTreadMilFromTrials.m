function CombFilesOfTreadMilFromTrials()
    filesLocationTxt  = '\\jackie-analysis\F\LayerV\Videos\PN11\02.03.22_ETL_Treadmill\t*.txt';
    outputpath = '\\jackie-analysis\F\LayerV\Videos\PN11\02.03.22_ETL_Treadmill\';
    doFixTwoP = false;
    BehaveFrameRate = 200;
    ImagingFrameRate = 30;
    timeT = 12;
    
    tfileTxt = dir(filesLocationTxt);
    
    treadmilDataAll = readtable(fullfile(tfileTxt(1).folder, tfileTxt(1).name));
    
    for i = 2:length(tfileTxt)
        temp  = readtable(fullfile(tfileTxt(i).folder, tfileTxt(i).name));
        diffTreadMil = temp.treadmill(1) -  treadmilDataAll.treadmill(end);
        diffTreadMilTime = temp.time(1) -  treadmilDataAll.time(end);
        diffTP = temp.twoP(1) -  treadmilDataAll.twoP(end);
        
        temp.treadmill = temp.treadmill - diffTreadMil;
        temp.time = temp.time - diffTreadMilTime; 
        temp.twoP = temp.twoP - diffTP;
        treadmilDataAll(end+1:(end + size(temp, 1)), :) = temp; 
    end
    
    if doFixTwoP
        countP = 0;
        for i = 1:(BehaveFrameRate / ImagingFrameRate):timeT*BehaveFrameRate*length(tfileTxt)
            treadmilDataAll.twoP(round(i):(round(BehaveFrameRate / ImagingFrameRate + i - 1))) = countP;
            countP = countP + 1;
        end
    end
    
    writetable(treadmilDataAll,[outputpath, '\AllFinal.txt']);
end