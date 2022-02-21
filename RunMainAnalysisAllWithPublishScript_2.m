addpath('C:\Users\Jackie\Desktop\Hadas\AnalysisProject\MatlabAnalysis');

globalParameters.TPAFolder = 'G:\Layer V\Analysis\';
globalParameters.MainFolder = '\\jackie-analysis\e\';
globalParameters.AnimalName = 'PT6';
globalParameters.DateAnimal = '07.15.21_Tuft_Control';
globalParameters.swcFile = 'swcFiles\neuron_2.swc';
globalParameters.neuronNumberName = 'N2';
globalParameters.RunnerDate = 'final';
globalParameters.RunnerNumber = 'Run1';
globalParameters.behaveType = 'no_behave';
globalParameters.analysisType = 'Pearson\SP';
globalParameters.treadmilFile = '';
globalParameters.trajFolderName = '';
globalParameters.videoPath = '';

globalParameters.doBehave = false;

globalParameters.swcLocationWithInOutputFolder = false;

globalParameters.neuronActivityPathTPA = fullfile(globalParameters.TPAFolder, globalParameters.AnimalName, globalParameters.DateAnimal);

globalParameters.neuronTreePathSWC = fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName, globalParameters.DateAnimal, globalParameters.swcFile);

globalParameters.behaveFileTreadMillPath = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, globalParameters.DateAnimal, globalParameters.treadmilFile);

globalParameters.eventsDetectionFolder = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, ...
    globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
    globalParameters.RunnerDate,globalParameters.RunnerNumber, 'EventsDetection');
mkdir(globalParameters.eventsDetectionFolder);

globalParameters.behaveTreadMilOutputFolder = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, ...
    globalParameters.DateAnimal, 'Analysis', 'BehaveTreadMilOutput');
    
globalParameters.isSimData = false;
globalParameters.reverseHeatMap = false;
globalParameters.hyperbolicDistMatrixLocation = "";
% globalParameters.hyperbolicDistMatrixLocation = fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName,...
%     globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
%     globalParameters.RunnerDate, globalParameters.RunnerNumber, 'HS_create\StructureTreeHS\matlab_matrixbernoulli_100_3000.mat'); 

globalParameters.glmResults = fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName,...
    globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'glmAnalysis_Running', 'glmResultsSummary.mat');

% -------------------------------------------------
globalParameters.runMantelTestPerDistanceThreshold = true;
globalParameters.runMantelTestPerDistanceThreshold_only = false;
globalParameters.MantelTJump = 50;
globalParameters.MantelTJumpFolder = fullfile(globalParameters.MainFolder, 'Shay',...
    'MantelThresholdTestSummary', globalParameters.AnimalName, globalParameters.DateAnimal,...
    globalParameters.neuronNumberName, globalParameters.RunnerDate,globalParameters.RunnerNumber, globalParameters.behaveType, globalParameters.analysisType);

% ----------------------------------------------------

%     Can be Euclidean OR ShortestPath OR ShortestPathCost OR Branch OR HyperbolicDist_L OR HyperbolicDist_P  = ( Between 2 roi according to the tree path )
globalParameters.roiTreeDistanceFunction = 'ShortestPath';
globalParameters.costSP = 0;

globalParameters.DistType = 'Pearson';

globalParameters.costCluster = -1;

globalParameters.std_treadMilThreshold = 0.0001;
globalParameters.aftertonetime = 26;

globalParameters.winLength = 2; % defaulte is 4;

globalParameters.splinesL = '\\jackie-analysis10\users\Jackie\Desktop\Hadas\AnalysisProject\MatlabAnalysis';

globalParameters.doExtraAnalysis = true;

globalParameters.treadMilExtraPlot = true;

globalParameters.reRunClusterData = false;
   
globalParameters.behaveFrameRateTM = 100;

globalParameters.ImageSamplingRate = 30;
globalParameters.time_sec_of_trial = 12;
globalParameters.trialNumber = [1, 2];
globalParameters.BehavioralSamplingRate = 200;
globalParameters.behavioralDelay = 0;
globalParameters.toneTime = 4;

%     No events put non
globalParameters.runByEvent = {'non'};
globalParameters.isHandreach = true;
globalParameters.EventTiming = 'start';
globalParameters.runByNegEvent = false;
globalParameters.FirstEventAfter = {'non'};

%     FOR Hand Reach 
%     NO labels 0, 1 suc , 2 fail
globalParameters.split_trialsLabel = 0;
globalParameters.runBehaveLag = [-inf, inf];
globalParameters.do_eventsBetween = false;
globalParameters.doBehaveAlignedPlot = false;

globalParameters.excludeTrailsByEventCount.Name = 'non';
globalParameters.excludeTrailsByEventCount.countRange = [-inf, inf];

%     Can be WindowEventFULLPearson OR
%     WindoEventToPeakPearson OR PeaksPearson OR glm
%     OR WindoEventToPeakCov OR WindoEventPearsonPerEvent OR HS_Activity
globalParameters.roiActivityDistanceFunction = 'WindoEventToPeakPearson';

globalParameters.clusterCount = 4;
globalParameters.eventWin = 10;

% 0.01 
globalParameters.mean_aV = 0.01; 
globalParameters.aVForAll = 0.3;
globalParameters.aVFix.location = [];
globalParameters.aVFix.values = [];

globalParameters.sigmaChangeValue = 0;
    
globalParameters.sigmaFix.location = [];
globalParameters.sigmaFix.values = [];

% most of the time 1.5 is good
globalParameters.thresholdGnValue = 1.5;
globalParameters.thresholdGnFix.location = [];
globalParameters.thresholdGnFix.values = [];

globalParameters.runMLS = true;

globalParameters.excludeRoi = [];

globalParameters.apical_roi = [];

globalParameters.firstDepthCompare = 1;
globalParameters.secDepthCompare = -1;

globalParameters.outputpath = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, ...
    globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
    globalParameters.RunnerDate,globalParameters.RunnerNumber, globalParameters.behaveType, globalParameters.analysisType);

globalParameters.outputpath = char(globalParameters.outputpath);
mkdir(globalParameters.outputpath);

globalParameters.mantelRF = fullfile(globalParameters.MainFolder, 'Shay', 'MantelSummary', globalParameters.AnimalName, globalParameters.DateAnimal,...
        globalParameters.neuronNumberName, globalParameters.RunnerDate, globalParameters.RunnerNumber, globalParameters.behaveType, globalParameters.analysisType);
mkdir(globalParameters.mantelRF);
    

warning('off','all');

codePar = 'PublishPDFMainAnalysis2(globalParameters)';

% publish('PublishPDFMainAnalysis2.m', 'showCode', false,...
%     'format','pdf', 'catchError', false, 'outputDir', globalParameters.outputpath, 'codeToEvaluate', codePar);

% PublishPDFMainAnalysis2(globalParameters);

% RunnerHadasCode(globalParameters);

mainRunnerNeuronTreeAndActivityAnalysis_V3(globalParameters);
    
clear;