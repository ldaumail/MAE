function LD_phantom_grating_compare_v12(subject, session, debug)

%In this version, we add multiple velocities
% subject = 'Dave';                                                                                                                                                                                                                                                     
% session = 1;                                                                                                                           
% debug = 1;


ex.version = 'v12';
% global EyeData rect w xc yc eye_used 
%%%% resolution 
if debug == 1
    % eyetracking on (1) or off (0)
    ET = 0;
    ex.screenWidth = 40;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 46;             % in cm; 3Tb/office=43, miniHelm=57;
	ex.resolution = SetResolution(max(Screen('Screens')),1280,1024,85); % laptop 1920,1080/ 2880, 1800 ,0
    ex.gammaCorrection = 0;       % make sure this = 1 when you're at the scanner!
else
    
    ET = 1;                                                                                                                             
    ex.screenWidth = 40;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 46;             % in cm; 3Tb/office=43, miniHelm=57;
    ex.resolution = SetResolution(max(Screen('Screens')),1280,1024,85); % ET room 1600,900,60
    ex.gammaCorrection = 1;       % make sure this = 1 when you're at the scanner!
end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
responseKeys = zeros(1,256);
responseKeys(KbName('0'))=1; % button box 1
responseKeys(KbName('1'))=1; % button box 1
responseKeys(KbName('2'))=1; % button box 2
responseKeys(KbName('3'))=1; % button box 1
responseKeys(KbName('4'))=1; % button box 1
responseKeys(KbName('5'))=1; % button box 1
responseKeys(KbName('Return'))=1; % button box 3
responseKeys(KbName('ENTER'))=1; % button box 3

Screen('Preference', 'SkipSyncTests', 0);

% ex.scanNum = input('Scan number :');
ex.runNum = input('Run number :');

%%% basic naming set-up
ex.subject = subject;
ex.session = session;


%%%% set-up rand
ex.startTime = clock;
rng(sum(100*ex.startTime));
ex.rand = rng;

%%%% files and things
ex.root = pwd;
ex.date = datestr(now,30);

%%%% 2D sine wave grating inducers properties
ex.stim.spatialFreqDeg = 0.5/2;   % cycles per degree of visual angle
ex.stim.orientation = [180]; %[90 180];                                                % in degrees
ex.stim.gaborHDeg = 12;                                                   % in degrees of visual angle
ex.stim.gaborWDeg = 6; 
ex.stim.gapSizeDeg = 4;
ex.stim.distFromFixDeg = (ex.stim.gapSizeDeg+ex.stim.gaborWDeg)/2;  %each grating center deg horizontal away from fixation (grating center 6 deg away)

ex.stim.backgroundLum = [60 60 60];
ex.stim.contrast = [0.03 0.03 0.15 0.15 0.60 0.60];
ex.stim.contrastOffset = [(ex.stim.backgroundLum(1)./255)./(1-ex.stim.contrast(1)), ex.stim.backgroundLum(1)./255, (ex.stim.backgroundLum(1,1)./255)./(1-ex.stim.contrast(3)),...
    ex.stim.backgroundLum(1)./255,(ex.stim.backgroundLum(1,1)./255)./(1-ex.stim.contrast(5)), ex.stim.backgroundLum(1)./255];%+ex.stim.contrast/2;
ex.stim.luminanceRange = 2*ex.stim.contrast.*ex.stim.contrastOffset;
ex.stim.contrastMultiplicator = ex.stim.luminanceRange./2;  % for sine wave

ex.stim.maxLum = 255*(ex.stim.contrastOffset+ex.stim.contrastMultiplicator);
ex.stim.minLum = 255*(ex.stim.contrastOffset-ex.stim.contrastMultiplicator);
ex.stim.contrast = (ex.stim.maxLum-ex.stim.minLum)./(ex.stim.maxLum+ex.stim.minLum);

%%%% sine wave grating timing (within block scale)
ex.initialFixation = 1;        % in seconds
ex.finalFixation = 2;          % in seconds
ex.blockLength = 10; %ex.trialFixation+ ceil(ex.stimDur*ex.stimsPerBlock);           % in seconds
ex.betweenBlocks = 2;          % in seconds

ex.flipsPerSec = 60;  % 60;         % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
ex.flipWin = 1/ex.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 
ex.stim.cycPerSec = 2; %drifting speed in cycles of grating per sec
ex.stim.motionRate = 360*ex.stim.cycPerSec; %drifting speed in degrees of visual angle per sec
ex.stim.dphase = ex.stim.motionRate/ex.flipsPerSec; %degrees per flip

% Reference grating stimulus
ex.test.spatialFreqDeg = 0.5/2; %0.286;                                          % cycles per degree of visual angle
ex.test.orientation = [180]; %[90 180];                                                % in degrees
ex.test.gaborHDeg = 12;                                                  % in degrees of visual angle
ex.test.gaborWDeg = 6;                                
ex.test.contrast = 0.03;
ex.test.contrastOffset = [(ex.stim.backgroundLum(1,1)./255)./(1+ex.test.contrast)]; %(ex.stim.backgroundLum(2,1)./255)./(1+ex.stim.contrast), ex.stim.backgroundLum(3,1)./255];%+ex.stim.contrast/2;
ex.test.luminanceRange = 2*ex.test.contrast*ex.test.contrastOffset;
ex.test.contrastMultiplicator = ex.test.luminanceRange./2;  % for sine wave


%%%% conditions & layout (across blocks scale)

ex.conds = {'LowContPhUp','LowContPhDown','LowContPhCtUp','LowContPhCtDown',...
    'MedContPhUp','MedContPhDown','MedContPhCtUp','MedContPhCtDown', ...%
    'HighContPhUp','HighContPhDown','HighContPhCtUp','HighContPhCtDown'
       }; 
ex.numConds = length(ex.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block
ex.repsPerRun = 4;              % repetitions of each condition per run
ex.nTrials = ex.numConds*ex.repsPerRun;
ex.numBlocks = ex.numConds*ex.repsPerRun;

ex.condShuffle = [];
ex.locShuffle = [];
for i =1:ex.repsPerRun
    ex.condShuffle = [ex.condShuffle, Shuffle([1:ex.numConds])];
    ex.locShuffle = [ex.locShuffle, Shuffle([1,2,1,2,1,2,1,2,1,2,1,2])]; %2 locations (left and right)
end
ex.locShuffle = repmat(ex.locShuffle,1,2);
%%%% fixation 
ex.fixSizeDeg =  .25;            % in degrees, the size of the biggest white dot in the fixation

%%%% screen
ex.backgroundColor = [127 127 127];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
ex.fontSize = 20;

%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

% %% create the timing model of stimulus conditions for this particular run
% clear i
for i =1:ex.numConds
    conditions(i).name = ex.conds(i);
    conditions(i).startTimes = [];
end

%%
%%%%%%%%%%%%%%%
% open screen %
%%%%%%%%%%%%%%%
HideCursor;
Priority(9);
%Priority(0);
%%%% open screen
screen=max(Screen('Screens'));
if debug
    [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat);
    %[w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[100 100 600 400],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
    
else
    %[w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[100 100 600 400],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
    
    [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
end
Screen(w, 'TextSize', ex.fontSize);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%gamma correction, file prepared for room 417
if ex.gammaCorrection
    % Gamma correction (run phase2_photometry.mat in 417C computer, get gamma
    % table)
    load("phase2_photometry.mat");
    Screen('LoadNormalizedGammaTable', screen, inverseCLUT);
end
%%%% timing optimization
frameInt = Screen('GetFlipInterval',w);
slack = frameInt/2;
frameRate =  1/frameInt;%Screen('NominalFrameRate',w);

flipTimes = [0:frameInt*frameRate/ex.flipsPerSec:ex.blockLength(1)];
flipTimes = flipTimes(1:length(flipTimes)-1);
ex.stim.flipTimes = flipTimes;

nconds = length(ex.stim.contrastMultiplicator);
ex.stim.phases = nan(nconds, length(ex.stim.flipTimes));
clear r l
for l =1:nconds
    ex.stim.phases(l,:) = ((1:length(ex.stim.flipTimes))-1)*ex.stim.dphase;%
end

%reasonning behind the calculation of ex.stim.phases:
%GOAL: render the back and forth of grating drifts oscillatory in time instead
%of linear to smooth the signal phase shifts at the time of drift direction
%reversal
%1) The length (x) of grating shift after each flip is the variable that will have to be oscillatory
%over time

%2)oscillatory signal formula is x(t) = A*cos(w*t+phase) (x is the oscillatory signal, w is the angular
%frequency in rad/s, t is the time in s, phase is a phase shift in rad at time t (here 0),
%A is the amplitude of the signal

%3)Here the angular frequency is calculated based on the idea of one
%oscillation per stimulus. Considering a stimulus as 1 cycle forward and 1 cycle backwards(definition might change based on how many cycles are desired per stimulus), one oscillation will represent 2 cycles =ex.stim.cycPerSec*ex.stimDur)
% thus the temporal frequency of the oscillation should be 1/ex.stimDur
% (this is w)

%%4) We multiply w by t. Here the resolution of the signal is given by the number of flips per second, so the duration is given by flipTimes.
%Finally, we multiply by 2*pi for the units to be in radian.

%5)Then we subtract the phase at each flip at time t which here is 0.

%6)The amplitude of cos() varies between -1 and +1. Here we want our
%displacement (grating shift) in degrees with a displacement of 360? (of the spatial grating) every half temporal oscillation so we multiply by 180 for the range to be [-180;+180].

%7) we multiply by the number of cycles desired to drift over one lap
%(ex.stim.cycles).

%%%% scale the stim params for the screen
ex.ppd = pi* rect(3) / (atan(ex.screenWidth/ex.viewingDist/2)) / 360;
ex.fixSize = round(ex.fixSizeDeg*ex.ppd);
ex.gapSize = round(ex.stim.gapSizeDeg*ex.ppd); %to draw the reference grating
ex.stim.distFromFix = round(ex.stim.distFromFixDeg*ex.ppd); %for the actual gap size
ex.gaborHeight = round(ex.stim.gaborHDeg*ex.ppd);                 % in pixels, the size of our objects
ex.gaborWidth = round(ex.stim.gaborWDeg*ex.ppd);                 % in pixels, the size of our objects 
ex.rawGaborHeight = ex.gaborHeight;
ex.rawGaborWidth = ex.gaborWidth;

%% Create a sinewave grating image saved for each drift phase and each condition

ex.rectSWave = nan(ex.rawGaborHeight,ex.rawGaborWidth,length(ex.stim.flipTimes), nconds);
ex.rectSWaveID = nan(length(ex.stim.flipTimes),nconds);
clear c r f
for c =1:nconds %
    for f = 1:length(ex.stim.flipTimes)
        phase = ex.stim.phases(c,f);
        ex.rectSWave(:,:,f,c) = makeSineGrating(ex.rawGaborHeight,ex.rawGaborWidth,ex.stim.spatialFreqDeg,...
            ex.stim.orientation,phase,ex.stim.contrastOffset(c),ex.stim.contrastMultiplicator(c),...
            ex.ppd);
        ex.rectSWaveID(f,c) = Screen('MakeTexture', w, squeeze(ex.rectSWave(:,:,f,c)));
    end
end

%% Creat low contrast sinewave

ex.testSWave = nan(ex.rawGaborHeight,ex.rawGaborWidth,length(ex.stim.flipTimes), nconds);
ex.testSWaveID = nan(length(ex.stim.flipTimes),nconds);
clear c f
for c =1:nconds %
    for f = 1:length(ex.stim.flipTimes)
    phase = ex.stim.phases(c,f)+180;
    ex.testSWave(:,:,f,c) = makeSineGrating(ex.rawGaborHeight,ex.rawGaborWidth,ex.stim.spatialFreqDeg,...
        ex.stim.orientation,phase,ex.test.contrastOffset,ex.test.contrastMultiplicator,...
        ex.ppd);
    ex.testSWaveID(f,c) = Screen('MakeTexture', w, squeeze(ex.testSWave(:,:,f,c)));
     end
end

%% Sine wave gratings locations (in the task loop since it changes)
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+e.vertOffset;

%% stim location
 xL = rect(3)/2-rect(3)/4; % % = stimulus located on the left side of the screen
 xR = rect(3)/2+rect(3)/4; % = stimulus located on the right side of the screen
% 
% yT = rect(4)/2 - (ex.stim.gapSizeDeg+ex.stim.gaborHDeg)*ex.ppd/2; % stimulus located 4 degrees above screen center
% yB = rect(4)/2+ (ex.stim.gapSizeDeg+ex.stim.gaborHDeg)*ex.ppd/2; % stimulus located 4 degrees below screen center


%left grating
xLL = xc-rect(3)/4 - ex.stim.distFromFix;
xLR = xc-rect(3)/4 + ex.stim.distFromFix;
% grating y locations
yL = yc;
%right grating
xRL = xc + rect(3)/4 - ex.stim.distFromFix;
xRR = xc + rect(3)/4 + ex.stim.distFromFix;
yR = yc;


% yC = rect(4)/2; % stimulus located on screen center
%% %%%% initial window - wait for backtick
DrawFormattedText(w,'Look at the blank gap and use the low contrast grating as an anchor  (feel free to look back and forth). Press: \n\n 0. If you don t see any pattern, \n\n 1. If you experience a very faint impression of grating pattern (does not connect all the way through) \n\n 2. Moderate impression of grating pattern (1/2 as strong as the physical grating) \n\n 3. Strong impression of grating pattern, but somewhat weaker than the physical grating \n\n 4. Vivid impression of grating pattern, as strong as the physical grating \n\n 5. Impression is stronger than the grating pattern \n\n Press Space to start'... % :  '...
    ,'center', 'center',[0 0 0]);
Screen(w, 'Flip', 0);
%WaitSecs(2);
KbTriggerWait(KbName('Space'), deviceNumber);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% START task TASK/FLIPPING
KbQueueCreate(deviceNumber,responseKeys);
n = 1; %drift phase
c = 1; %condition
% t =1; %trial number
blockCnt = 1;
cntCond = zeros(length(ex.conds),1);
ex.responseTimes=[];
ex.resp = [];
ex.trialStartTime = [];
ex.condResp = [];
%%% initial fixation
if n == 1 && blockCnt == 1 %for first block
    ex.tasktstart = clock;
    ex.startRun = GetSecs();
    Screen('FillRect', w, ex.stim.backgroundLum);
    Screen(w, 'Flip', 0);
    WaitSecs(ex.initialFixation);
end
%%% Launch the task
while(1) %n <= length(ex.trialFlips)
    KbQueueStart();
    condNum = ex.condShuffle(c);
    condName = conditions(condNum).name{:}; 
    l = ceil(condNum/2); % indices for each contrast level

    Screen('FillRect', w, ex.stim.backgroundLum)
    %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ex.locShuffle(c) == 1
        ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xLL,yL);
        ex.rectRRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xLR,yL);
        ex.rectCRect =  CenterRectOnPoint([0 0 ex.gapSize ex.rawGaborHeight],xR,yc);
    elseif ex.locShuffle(c) == 2
        ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xRL,yR);
        ex.rectRRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xRR,yR);
        ex.rectCRect =  CenterRectOnPoint([0 0 ex.gapSize ex.rawGaborHeight],xL,yc);
    end
    % stim
%     Screen('DrawTexture', w, ex.rectSWaveID(n,thisCond),[],ex.rectLRect);
%     Screen('DrawTexture', w, ex.rectSWaveID(n,thisCond),[],ex.rectRRect);
%     
    % stim
    if contains(condName, 'Up')
        Screen('DrawTexture', w, ex.rectSWaveID(n,l),[],ex.rectLRect);
        Screen('DrawTexture', w, ex.rectSWaveID(n,l),[],ex.rectRRect);
        Screen('DrawTexture', w, ex.testSWaveID(n,l),[],ex.rectCRect);
        
    elseif contains(condName, 'Down')
        Screen('DrawTexture', w, ex.rectSWaveID(end-(n-1),l),[],ex.rectLRect);
        Screen('DrawTexture', w, ex.rectSWaveID(end-(n-1),l),[],ex.rectRRect);
        Screen('DrawTexture', w, ex.testSWaveID(end-(n-1),l),[],ex.rectCRect);
        
    end

    
    DrawFormattedText(w,' 0. No pattern, 1. Very faint impression 2. Moderate impression (50%) 3. Strong impression  4. As strong as the physical grating \n\n 5. Stronger than the pattern'... % :  '...
    ,20,970,[0 0 0]);
    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n == 1
        [VBLT, ex.startTrial, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        flipTimes = ex.startTrial;
        ex.trialStartTime = [ex.trialStartTime, ex.startTrial];
        ex.condResp = [ex.condResp, c];
    else
        [VBLT,flipTime, FlipT, missed] = Screen(w, 'Flip',ex.startTrial + ex.stim.flipTimes(n) - slack); %,   %%% ex.flipTime(n,c)
        flipTimes = [flipTimes, flipTime];
        
    end
    
    
    KbQueueStop();
    [pressed, firstPress]= KbQueueCheck();
    if (pressed == 1) && (firstPress(KbName('0')) > 0 || firstPress( KbName('1')) > 0 || firstPress(KbName('2')) > 0 || firstPress(KbName('3')) > 0 || firstPress(KbName('4')) > 0 || firstPress(KbName('5')) > 0) %%

            if (firstPress(KbName('0')) > 0)
                ex.resp = [ex.resp, 0];
                ex.responseTimes = [ex.responseTimes, firstPress(KbName('0')) - ex.startRun];
                
            elseif (firstPress(KbName('1')) > 0)
                ex.resp = [ex.resp, 1];
                ex.responseTimes = [ex.responseTimes, firstPress(KbName('1')) - ex.startRun];   
            elseif (firstPress(KbName('2')) > 0)
                ex.resp = [ex.resp, 2];
                ex.responseTimes = [ex.responseTimes, firstPress(KbName('2')) - ex.startRun];
            elseif (firstPress(KbName('3')) > 0)
                ex.resp = [ex.resp, 3];
                ex.responseTimes = [ex.responseTimes, firstPress(KbName('3')) - ex.startRun];
            elseif (firstPress(KbName('4')) > 0)
                ex.resp = [ex.resp, 4];
                ex.responseTimes = [ex.responseTimes, firstPress(KbName('4')) - ex.startRun];
            elseif (firstPress(KbName('5')) > 0)
                ex.resp = [ex.resp, 5];
                ex.responseTimes = [ex.responseTimes, firstPress(KbName('5')) - ex.startRun];
            end
            pressed = 0;

    end
     KbQueueFlush();
    n = n+1;
    if (n == length(ex.rectSWaveID(:,1,1))+1) % for any other block , reset frame index when previous trial ends
        n = 1;
        
    elseif (pressed && ismember(find(firstPress,1), [KbName('Return') KbName('ENTER')]))
        n = 1;
        c = c+1;
    end
    if c > ex.nTrials
        break;
    end
end



%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

ex.runTime = GetSecs - ex.startRun;
ex.rectSWave = [];
ex.testSWave = [];
savedir = fullfile(ex.root,'data',sprintf('dyn_MAE/s%s_%s/',subject,ex.version));
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(sprintf('/s%s_percept_test_%s_date%s_fix',subject,ex.version,num2str(ex.date)), '.mat'));
%save(savename,'ex');
save(savename,'ex','-v7.3')

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;                                                                                                                          


end 