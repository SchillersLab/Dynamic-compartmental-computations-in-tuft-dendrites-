function createVideoForTreadMil
    activityDataLocation = '\\jackie-analysis\e\Shay\4481\2019-03-07_4481_motor_aligned\Analysis\N3\Structural_VS_Functional\final\Run2\no_behave\Pearson\SP\roiActivityRawData.mat';
    behaveVideoLocation = 'F:\2019-03-07_4481\2019-03-07_4481_motor_behavior\trial0001.mkv';
    outputpath = '\\jackie-analysis\e\Shay\4481\2019-03-07_4481_motor_aligned\Analysis\N3\';
    treadmilTxtFile = 'F:\2019-03-07_4481\2019-03-07_4481_motor_behavior\trial0001_treadmill.txt';
    behaveTreadMil = '\\jackie-analysis\e\Shay\4481\2019-03-07_4481_motor_aligned\Analysis\BehaveTreadMilOutput\BehaveTreadMilResults.mat';
    treadmilData = readtable(treadmilTxtFile);
   
    mkdir(outputpath);
    
    v_behave = VideoReader(behaveVideoLocation);
    load(activityDataLocation, 'roiActivity');
    meanRoiActivity = mean(roiActivity, 2);
    
    load(behaveTreadMil, 'BehaveDataTreadmil');
    speedB = BehaveDataTreadmil.speed;
    
    totalFrames = length(BehaveDataTreadmil.walkconstant) + length(BehaveDataTreadmil.walkacceleration);
    
    vOut = VideoWriter([outputpath, '\behaveWithActivity2_2.avi']);    
    
    frameIndex = sort([BehaveDataTreadmil.walkconstant; BehaveDataTreadmil.walkacceleration]);
    open(vOut);

%     videoPlayer = vision.VideoPlayer;
        
    k = 1;
    lastLocationInTwoP = 0;
    treadIndex = 1;
    
    while(hasFrame(v_behave))
        if k > length(meanRoiActivity)
            break;
        end
        
        if any(frameIndex == k)
            fig = figure('visible','off'); hold on;
            plot(meanRoiActivity + 1, 'Color', 'k', 'LineWidth', 1.5);
            plot(speedB(1:size(meanRoiActivity, 1)) + 2, 'Color', [25 110 180, 255] ./ 255);

            fig.Position = [fig.Position(1), fig.Position(2), v_behave.Width, v_behave.Height];

            scatter(k, meanRoiActivity(k)+ 1, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'SizeData', 15);
            scatter(k, speedB(k) + 2, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'SizeData', 15);
            xlim([max(k - 100, 0), min(k +100, size(meanRoiActivity, 1))]);
            ylim([-1, 8]);
            % Save the frame in structure for later saving to video file

            cur = readFrame(v_behave);
            fin = getframe(gcf);

            imgt = horzcat(fin.cdata, cur);

    %         step(videoPlayer, imgt);
            writeVideo(vOut,imgt)
            
            close gcf
            clear fig imgt fin cur;
        else
            readFrame(v_behave);
        end
        
        if (treadmilData.twoP(treadIndex) ~= lastLocationInTwoP)
            k = k+1;
            lastLocationInTwoP = treadmilData.twoP(treadIndex);
        end
        
        treadIndex = treadIndex + 1;     
    end

%     release(videoPlayer);
    close(vOut)
    close all
end