function RunnerHadasCode(globalParameters)
    addpath('\\jackie-analysis10\users\Jackie\Desktop\Hadas\AnalysisProject\MatlabAnalysis');

    xmlfile = fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName, globalParameters.DateAnimal, 'XmlByBoth.xml');
    
    folderAnimal = fullfile(globalParameters.TPAFolder, globalParameters.AnimalName);
    trajpth = fullfile('\\192.114.21.82\g\Layer V\Videos', globalParameters.AnimalName , globalParameters.trajFolderName, 'DLC');
    roiListNamesPath = fullfile(globalParameters.outputpath, 'roiActivityRawData.mat');
    
    folderAnimalOutputPath= fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName);
    
    specific_experiment = globalParameters.DateAnimal;

    listExperiments = dir (folderAnimal);

    for index = 1: length(listExperiments)      
        if ~strcmp(listExperiments(index).name, '..') && ~strcmp(listExperiments(index).name, '.')

            if strcmp(specific_experiment, '') || (~strcmp(specific_experiment, '') && strcmp(listExperiments(index).name, specific_experiment))

                bda_tpa_folder = strcat(listExperiments(index).folder, '\',listExperiments(index).name);
                listFiles = dir(bda_tpa_folder);

                if ~isempty(listFiles)
                    bdaCount = 1;
                    for i = 1: length(listFiles)
                        testBDA = listFiles(i).name;
                        if contains({testBDA}, 'BDA')
                            BdaTpaList(bdaCount).BDA = [bda_tpa_folder '\' testBDA]; 

                            for k = 1: length(listFiles)
                                if contains({listFiles(k).name}, 'TPA')
                                    testTPA = strrep(listFiles(k).name,'TPA','BDA');
                                    if (strcmp(testTPA, testBDA))
                                        BdaTpaList(bdaCount).TPA = [bda_tpa_folder '\' listFiles(k).name]; 
                                        bdaCount = bdaCount + 1;
                                    end
                                end
                            end
                        end
                    end
                    BdaTpaList = getTrajFiles(BdaTpaList, trajpth, listExperiments(index).name);

                    BdaTpaList(1).roiListNamesPath = roiListNamesPath;
                    BdaTpaList(1).predictor = [];
                    
                    outputPath = strcat(folderAnimalOutputPath , '\' , listExperiments(index).name ,'\Analysis\', globalParameters.neuronNumberName ,'\atogram_printer');
                    mkNewFolder(outputPath);
                    runAnalysis(outputPath, xmlfile, BdaTpaList, 'atogram_printer', 'RoiSplit');

                    outputPath = strcat(folderAnimalOutputPath , '\' , listExperiments(index).name ,'\Analysis\', globalParameters.neuronNumberName ,'\delay2events');
                    mkNewFolder(outputPath);
                    runAnalysis(outputPath, xmlfile, BdaTpaList, 'delay2events', 'RoiSplit');

                    outputPath = strcat(folderAnimalOutputPath , '\' , listExperiments(index).name ,'\Analysis\', globalParameters.neuronNumberName, '\pcaTrajectories');
                    mkNewFolder(outputPath);
                    runAnalysis(outputPath, xmlfile, BdaTpaList, 'pcaTrajectories', 'RoiCorrelation');

                    outputPath = strcat(folderAnimalOutputPath , '\' , listExperiments(index).name ,'\Analysis\', globalParameters.neuronNumberName, '\centralityAnalysis');
                    mkNewFolder(outputPath);
                    runAnalysis(outputPath, xmlfile, BdaTpaList, 'centralityAnalysis', 'RoiCorrelation');

                    if globalParameters.isHandreach
                        outputPath = strcat(folderAnimalOutputPath , '\' , listExperiments(index).name ,'\Analysis\', globalParameters.neuronNumberName ,'\glmAnalysis');
                        mkNewFolder(outputPath);
                        runAnalysis(outputPath, xmlfile, BdaTpaList, 'glmAnalysis', 'RoiSplit');
                    end
                    
                    close all;
                    BdaTpaList = [];
                end
            end
        end
    end
end