MainFolder = '\\jackie-analysis\e\';
AnimalName = 'SM06';
DateAnimal = '02.03.21_Tuft_Treadmill_FreeRun_Trials\';
neuronNumberName = 'N1';
RunnerDate = 'final';
RunnerNumber = 'Run1_MLS';
behaveType = 'no_behave';
analysisType = 'Pearson\SP';

globalParamsFileLocation = fullfile(MainFolder, 'Shay' , AnimalName, ...
    DateAnimal, 'Analysis', neuronNumberName, 'Structural_VS_Functional',...
    RunnerDate,RunnerNumber, behaveType, analysisType, 'runParametes.mat');

load(globalParamsFileLocation);

globalParameters.outputpath = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, ...
    globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
    globalParameters.RunnerDate,globalParameters.RunnerNumber, globalParameters.behaveType, 'glm\SP');

globalParameters.doExtraAnalysis = false;
globalParameters.treadMilExtraPlot = false;
globalParameters.analysisType = 'glm\SP';
globalParameters.DistType = 'glm';
globalParameters.roiActivityDistanceFunction = 'glm';

globalParameters.runMantelTestPerDistanceThreshold = false;
globalParameters.runMantelTestPerDistanceThreshold_only = false;
globalParameters.MantelTJump = 50;
globalParameters.MantelTJumpFolder = fullfile(globalParameters.MainFolder, 'Shay',...
'MantelThresholdTestSummary', globalParameters.AnimalName, globalParameters.DateAnimal,...
globalParameters.neuronNumberName);
globalParameters.sigmaChangeValue = 0;
globalParameters.sigmaFix.location = [];
globalParameters.sigmaFix.values = [];
% most of the time 1.5 is good
globalParameters.thresholdGnValue = 4;
globalParameters.thresholdGnFix.location = [];
globalParameters.thresholdGnFix.values = [];
globalParameters.runMLS = true;
globalParameters.reverseHeatMap = false;

mainRunnerNeuronTreeAndActivityAnalysis_V3(globalParameters);
clear;