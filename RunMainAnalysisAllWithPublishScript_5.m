addpath('\\jackie-analysis10\users\Jackie\Desktop\Hadas\AnalysisProject\MatlabAnalysis');

globalParameters.TPAFolder = '\\jackie-analysis\e\ShayCode\pythonProject\larkumEtAl2009_2\';
globalParameters.MainFolder = '\\jackie-analysis\e\';
globalParameters.AnimalName = 'Simulation\4481_N3\Cai_1308\';
globalParameters.DateAnimal = '05.01.SpecificInputs_Control_Sca0.6_time750ms_BacK80BackFR1_run8_0_0_SR50';

% globalParameters.swcFile = '\..\Swc\cell_N3_4481_allseg_onlyTuft_7-3d_noObliq_subP.swc';

% globalParameters.swcFile = '\..\Swc\cell_N3_4481smTuft2_new_onlyTuft_7-3d_noObliq_subP.swc';
globalParameters.swcFile = '\..\Swc\cell_N3_4481_onlyTuft_7-3d_noObliq_subP.swc';

% globalParameters.swcFile = '\..\Swc\cell_N6_SM01_onlyTuft_7-3d_noObliq_subP.swc';

% globalParameters.swcFile = '\..\Swc\cell_N1_SM05_onlyTuft_7-3d_noObliq_subP.swc';
% globalParameters.swcFile = '\..\Swc\cell_N1_CT14_onlyTuft_7-3d_noObliq_subP.swc';
%  globalParameters.swcFile = '\..\Swc\cell_N1_SM03_onlyTuft_6-3d_noObliq_subP.swc';
 
globalParameters.neuronNumberName = 'N1';
globalParameters.RunnerDate = 'final';
globalParameters.RunnerNumber = 'Run1';
globalParameters.behaveType = 'no_behave';                                                                                                                                                        veType = 'no_behave';
globalParameters.analysisType = 'Pearson\SP';
globalParameters.treadmilFile = '';
globalParameters.trajFolderName = '';
globalParameters.videoPath = '';

globalParameters.hyperbolicDistMatrixLocation = "";
% globalParameters.hyperbolicDistMatrixLocation = fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName,...
%     globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
%     globalParameters.RunnerDate, globalParameters.RunnerNumber, 'HS_create\StructureTreeHS\matlab_matrixbernoulli_100_3000.mat'); 

globalParameters.glmResults = fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName,...
    globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'glmAnalysis', 'glmResultsSummary.mat');

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
globalParameters.doZscoreAll = false;
globalParameters.costCluster = -1;

globalParameters.std_treadMilThreshold = 0.001;
globalParameters.aftertonetime = 26;

globalParameters.winLength = 2; % defaulte is 4;

globalParameters.splinesL = 'C:\Users\Jackie\Desktop\HadasCode\AnalysisProject\MatlabAnalysis';

globalParameters.doExtraAnalysis = false;

globalParameters.treadMilExtraPlot = false;

globalParameters.reRunClusterData = false;
   
globalParameters.behaveFrameRateTM = 100;

globalParameters.ImageSamplingRate = 50;
globalParameters.time_sec_of_trial = 0.75;
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
globalParameters.mean_aV = 0.05; %0.05 for reg 0.01, for simData 0.8, for V data 0.08
globalParameters.aVForAll = 0.3;
globalParameters.aVFix.location = [];
globalParameters.aVFix.values = [];

globalParameters.sigmaChangeValue = 0;
    
globalParameters.sigmaFix.location = [];
globalParameters.sigmaFix.values = [];

% most of the time 1.5 is good
globalParameters.thresholdGnValue = 0.05; %0.05 sim , V 3
globalParameters.thresholdGnFix.location = [];
globalParameters.thresholdGnFix.values = [];

globalParameters.runMLS = false;
globalParameters.isSimData = false;

% 4481
globalParameters.reverseHeatMap = true;

% SM03 & 4482
% globalParameters.reverseHeatMap = true;

% 4481_N3
% globalParameters.excludeRoi = [3457, 2347, 242, 925, 721, 4560, 3798, 4473, 1279, 879, 10002, 5572,4041,4923,4814,5127,3657,2546,2442,2652,2853,3067,1166,1272,1388,5342,1099,967,1439,1685,1776,2171,2041,533,295,290];
% 
% SM03_N1
% globalParameters.excludeRoi = [812,463,221,235,448,411,203,931,720,975,1089,1211,1252,1297,1424,1499,1557,1724,2844,1869,2539,2320,1879,1955,2121,2067,3052,3160,3277,3433,3474,3677,3692,3800];

% 4482_N1
% globalParameters.excludeRoi = [2786,3215,3061,3463,1503,1871,2604,1370,3537,2302,2706,1049,998,639,952,628,480,333];

% CT14|_N1
% globalParameters.excludeRoi = [5152,5788,6290,6043,5798,5988,5767,5208,5492,4960,4728,4837,4010,756,6378,4013,4631,4037,4164,1217,1621,1477,1069,470,6657,4493,5076,2961,3968,3512,3030,3363,2507,2588,2900,2737,2670,2725,2387,2083];

globalParameters.excludeRoi = [];

globalParameters.apical_roi = [];

globalParameters.firstDepthCompare = 1;
globalParameters.secDepthCompare = -1;

globalParameters.outputpath = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, ...
    globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
    globalParameters.RunnerDate,globalParameters.RunnerNumber, globalParameters.behaveType, globalParameters.analysisType);

globalParameters.outputpath = char(globalParameters.outputpath);
mkdir(globalParameters.outputpath);

warning('off','all');


% codePar = 'PublishPDFMainAnalysis3(globalParameters)';
% 
% publish('PublishPDFMainAnalysis5.m', 'showCode', false,...
%     'format','pdf', 'catchError', false, 'outputDir', globalParameters.outputpath, 'codeToEvaluate', codePar);
% 
% PublishPDFMainAnalysis4(globalParameters);

% RunnerHadasCode(globalParameters)

mainRunnerNeuronTreeAndActivityAnalysis_V3(globalParameters);

clear;
