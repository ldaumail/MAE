function LD_static_MAE_revisions(subject, session, debug)

%In this version, we add multiple velocities
% subject = 'sub-01'; 
% session = 1;                                                                                                                           
% debug = 1;


ex.version = 'v1';
%%%% resolution 
if debug == 1

    ex.screenWidth = 40;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 46;             % in cm; 3Tb/office=43, miniHelm=57;
	ex.resolution = SetResolution(max(Screen('Screens')),3024,1964,120); % laptop 1920,1080/ 2880, 1800 ,0
    ex.gammaCorrection = 0;       % make sure this = 1 when you're at the scanner!
else                                                                                                                            
    ex.screenWidth = 40;             % in cm; % 16 in eye tracking room 425%laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 46;             % in cm; %23 in eye tracking                                                                                                                          room 425 3Tb/office=43, miniHelm=57;
    ex.resolution = SetResolution(max(Screen('Screens')),1280,1024,85); % ET room 1600,900,60
    ex.gammaCorrection = 1;       % make sure this = 1 when you're at the scanner!
end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
responseKeys = zeros(1,256);
responseKeys(KbName('LeftArrow'))=1; % button box 
responseKeys(KbName('RightArrow'))=1; % button box 
responseKeys(KbName('DownArrow'))=1; % button box 3

responseKeys(KbName('Space'))=1; % button box 3
%Screen('Preference', 'SkipSyncTests', 0);
Screen('Preference', 'SkipSyncTests', 1);
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
ex.stim.spatialFreqDeg = 0.5/2; %0.36;%  % cycles per degree of visual angle
ex.stim.orientation = [90]; %[90 180];                                                % in degrees
ex.stim.gaborHDeg = 6; %5.2;                                                   % in degrees of visual angle
ex.stim.gaborWDeg = 12; %8.8; 
ex.stim.gapSizeDeg = 2.6; %4; 
ex.stim.distFromFixDeg = (ex.stim.gaborHDeg+ex.stim.gapSizeDeg)/2;%3;%2; %1.5 %each grating edge 1.5 deg horizontal away from fixation (grating center 6 deg away)

ex.stim.backgroundLum = [60 60 60];
ex.stim.contrast = [0.075 0.075 0.15 0.15 0.60 0.60];
ex.stim.contrastOffset = [(ex.stim.backgroundLum(1)./255)./(1-ex.stim.contrast(1)), (ex.stim.backgroundLum(1)./255)./(1-ex.stim.contrast(1)),(ex.stim.backgroundLum(1,1)./255)./(1-ex.stim.contrast(3)),...
    (ex.stim.backgroundLum(1,1)./255)./(1-ex.stim.contrast(3)),(ex.stim.backgroundLum(1,1)./255)./(1-ex.stim.contrast(5)), (ex.stim.backgroundLum(1,1)./255)./(1-ex.stim.contrast(5))];%+ex.stim.contrast/2;
ex.stim.luminanceRange = 2*ex.stim.contrast.*ex.stim.contrastOffset;
ex.stim.contrastMultiplicator = ex.stim.luminanceRange./2;  % for sine wave

ex.stim.maxLum = 255*(ex.stim.contrastOffset+ex.stim.contrastMultiplicator);
ex.stim.minLum = 255*(ex.stim.contrastOffset-ex.stim.contrastMultiplicator);
ex.stim.contrast = (ex.stim.maxLum-ex.stim.minLum)./(ex.stim.maxLum+ex.stim.minLum);



%%%% sine wave grating timing (within block scale)
ex.initialFixation = 6;        % in seconds
ex.finalFixation = 2;          % in seconds
ex.blockLength = 3; %45; %ex.trialFixation+ ceil(ex.stimDur*ex.stimsPerBlock);           % in seconds
ex.ITI2 = .5;% in seconds
% ex.betweenBlocks = 2;          % in seconds
ex.flipsPerSec = 60;  % 60;         % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
ex.flipWin = 1/ex.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 
ex.stim.cycPerSec = 1.5; % 
ex.stim.motionRate = 360*ex.stim.cycPerSec; %drifting speed in degrees of visual angle per sec
ex.stim.dphase = ex.stim.motionRate/ex.flipsPerSec; %degrees per flip


%%%% Test stimulus: counterphasing grating
ex.test.spatialFreqDeg = ex.stim.spatialFreqDeg;
ex.test.contrast = 0.075;
ex.test.contrastOffset = ex.stim.backgroundLum(1,1)./255; %(ex.stim.backgroundLum(1,1)./255)./(1-ex.test.contrast);%ex.stim.backgroundLum(:,1)./255;% 
ex.test.luminanceRange = 2*ex.test.contrast*ex.test.contrastOffset;%0.1; %linspace(0.01,0.20,10);%[0.05, 0.10, 0.15];                                                 % in %, maybe?? %here the number of stimulus contrast levels is the number of different conditions
ex.test.contrastMultiplicator = ex.test.luminanceRange/2;  % for sine wave 0.5 = 100% contrast, 0.2 = 40%

%0.2+0.7/2;%0.425;
ex.test.gaborHDeg = (2/3)*ex.stim.gapSizeDeg; %ex.stim.gapSizeDeg*(2/3);                                                  % in degrees of visual angle
ex.test.gaborWDeg = ex.stim.gaborWDeg; %ex.stim.gaborWDeg;
ex.test.distFromFixDeg = 0; % in degrees of visual angle, grating center 2 deg away (edge 1 deg away)
ex.test.cycPerSec = ex.stim.cycPerSec;
ex.testDur = (1/ex.test.cycPerSec);        % in seconds. 1.77 sec refers to sine wave grating 1.77 = 2cycles/1.13cyc.sec-1 mutiplied by 2 for back and forth

%%%% Fixation
% ex.fixSizeDeg =  .2;            % in degrees, the size of the biggest white dot in the fixation
% ex.bigFixSizeDeg = 0.5;
% ex.outerFixPixels = 2;          % in pixels, the black ring around fixation

%%%% Nonius lines 
ex.lineHdeg = 0.45;
ex.lineWdeg = 0.06;

%%% horizontal line 
ex.horiLineWdeg = 0.45;


%%%% global stim/test offset from fixation

ex.stim.horiDistFromFixDeg = ex.stim.gaborWDeg*3/8;

%%%% conditions & layout (across blocks scale)

ex.conds = {'LowContPhRightAlt','LowContPhLeftAlt', ...
    'LowContPhRightSim','LowContPhLeftSim', ...
    'MedContPhRightAlt','MedContPhLeftAlt', ...%
    'MedContPhRightSim','MedContPhLeftSim', ...%
    'HighContPhRightAlt','HighContPhLeftAlt', ...
    'HighContPhRightSim','HighContPhLeftSim'}; 
ex.repsPerRun = repmat(3,length(ex.conds),1);%3 trials per motion direction             % repetitions of each condition per run
condIdx = 1:length(ex.conds); %[1,4,7]; %conditions we are interested to keep
ex.conds = ex.conds(condIdx);
ex.repsPerRun = ex.repsPerRun(condIdx);
ex.numConds = length(ex.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block

ex.numBlocks = sum(ex.repsPerRun);
ex.condShuffle = [];
for i =1:ex.repsPerRun(1)
   ex.condShuffle = [ex.condShuffle, Shuffle(1:length(ex.conds))];
end


ex.trialFlips = (0:ex.flipWin:ex.blockLength(1));%+ex.testLength
ex.trialFlips = ex.trialFlips(1:end-1);

%%%% screen
ex.backgroundColor = [127.4933  127.4933  127.4933];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
ex.fontSize = 26;

%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

ex.onSecs = [ones(1,ex.blockLength(1))]; %zeros(1,ex.testLength)
ex.longFormBlocks = Expand(ex.onSecs,ex.flipsPerSec,1); %1 when block, 0 when between block
length(ex.longFormBlocks)
ex.stimOnSecs = [ones(1,ex.blockLength(1))]; % zeros(1,ex.testLength) 
ex.longFormStimOnSecs = Expand(ex.stimOnSecs,ex.flipsPerSec,1); %1 when stim on, 0 when fixation or between blocks

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
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % Set up alpha-blending for smooth (anti-aliased) lines

%gamma correction, file prepared for room 425
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
% 
flipTimes = [0:1/ex.flipsPerSec:1/(ex.stim.cycPerSec)]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames
ex.stim.flipTimes = flipTimes(1:length(flipTimes)-1);

nconds = length(ex.stim.contrastMultiplicator);
ex.stim.phases = nan(nconds, length(ex.stim.flipTimes));

clear r l

for l =1:nconds
    
    ex.stim.phases(l,:) = ((1:length(ex.stim.flipTimes))-1)*ex.stim.dphase;%(ex.stim.oscillation1(c,l,r,:).*180*ex.stim.cycles(1)+ ex.stim.oscillation2(c,l,r,:).*180*ex.stim.cycles(2))/2 + ex.stim.spatialPhase; %./ex.stimDur-2*pi*flipTimes./ex.stimDur make it oscillatory
    
end


flipTimesTest = [0:1/ex.flipsPerSec:1];%till 1 sec %multiply frameInt by 60/12 = 5 to flip the image every 5 frames
ex.test.flipTimes = flipTimesTest(1:length(flipTimesTest)-1);

ex.test.tempPhase1 = 0;% rand(1,1)*2*pi;
halfPeriod1 = linspace(-1,1,ex.testDur*ex.flipsPerSec/2+1);
halfPeriod2 = linspace(1,-1,ex.testDur*ex.flipsPerSec/2+1);
ex.test.oscillation1 = [halfPeriod1(1:end-1),halfPeriod2(1:end-1)];%cos(2*pi*(1/ex.testDur(1))*ex.test.flipTimes+ex.test.tempPhase1);
ex.test.spatialPhase1 = 0;
ex.test.tmpphase = ex.test.oscillation1(1).*90+ex.test.spatialPhase1; %./ex.stimDur-2*pi*flipTimes./ex.stimDur make it oscillatory

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

ex.gaborHeight = round(ex.stim.gaborHDeg*ex.ppd);                 % in pixels, the size of our objects
ex.gaborWidth = round(ex.stim.gaborWDeg*ex.ppd);                 % in pixels, the size of our objects
ex.rawGaborHeight = ex.gaborHeight;
ex.rawGaborWidth = ex.gaborWidth;
ex.stim.distFromFix = round(ex.stim.distFromFixDeg*ex.ppd);
ex.stim.gapSize = round(ex.stim.gapSizeDeg*ex.ppd);
ex.stim.fullHeight = round((ex.stim.gapSizeDeg+ex.rawGaborHeight*2)*ex.ppd);

%%%% scale the test params for the screen
ex.probeHeight = round(ex.test.gaborHDeg*ex.ppd);                 % in pixels, the size of our objects
ex.probeWidth = round(ex.test.gaborWDeg*ex.ppd);                 % in pixels, the size of our objects
ex.rawProbeHeight = ex.probeHeight*1;
ex.rawProbeWidth = ex.probeWidth*1;


%%%scale the fixation params for the screen
% ex.fixSize = round(ex.fixSizeDeg*ex.ppd);
% ex.bigFixSize = round(ex.bigFixSizeDeg*ex.ppd);
ex.lineW = round(ex.lineWdeg*ex.ppd);
ex.lineH = round(ex.lineHdeg*ex.ppd);
ex.horiLineW = round(ex.horiLineWdeg*ex.ppd);

%% Create sinewave grating frames saved for each repetition and each condition

ex.rectSWave = nan(ex.rawGaborHeight,ex.rawGaborWidth,length(ex.test.oscillation1),nconds);
ex.rectSWaveID = nan(length(ex.test.oscillation1),nconds);
clear c r

for l = 1:length(ex.stim.contrastMultiplicator)
    for f =1:length(ex.test.oscillation1)
        phase = ex.stim.phases(l,f);
        ex.rectSWave(:,:,f,l) = makeSineGrating(ex.rawGaborHeight,ex.rawGaborWidth,ex.stim.spatialFreqDeg,...
            ex.stim.orientation,phase,ex.stim.contrastOffset(l),ex.stim.contrastMultiplicator(l),...
            ex.ppd);
        ex.rectSWaveID(f,l) = Screen('MakeTexture', w, squeeze(ex.rectSWave(:,:,f,l)));
    end
end

% for i =1:60
% figure();
% imshow(squeeze(ex.rectSWave(:,:,i,1))/255)
% end
%check luminances ranges
% minval = min(squeeze(ex.rectSWave(1,1,1,:,:)),[],'all');
% maxval = max(squeeze(ex.rectSWave(1,1,1,:,:)),[],'all');

%% create dynamic grating image as a test for other eye
spphase = 0;
ex.testStim = makeCounterPhasingGrating(ex.rawProbeHeight,ex.rawProbeWidth,ex.test.spatialFreqDeg,...
    ex.stim.orientation,spphase,ex.test.contrastOffset(1),ex.test.contrastMultiplicator(1), ex.test.tmpphase(1),...
    ex.ppd);
ex.testStimID = Screen('MakeTexture', w, squeeze(ex.testStim));

%% Create random lines background for enhancing disparity effect
ex.bg.ih = ex.rawGaborHeight*2+2*ex.stim.gapSize;
ex.bg.iw = ex.rawProbeWidth*5/4;
ex.bg.freq = 11; % line segments per square degree
ex.bg.angle = 45;
ex.bg.len = 16; %segment length in number of pixels 
ex.background = makeRandomLines(ex.bg.ih,ex.bg.iw,ex.bg.freq,ex.bg.angle,ex.bg.len,ex.viewingDist,ex.stim.backgroundLum(1,:), ex.ppd);
ex.backgroundID = Screen('MakeTexture', w, ex.background);

%% Create blank gap intervening background
ex.blankGap = ones(ex.rawGaborHeight, ex.stim.gapSize).*ex.stim.backgroundLum(1,1);
ex.blankGapID = Screen('MakeTexture', w, ex.blankGap);

%% Create test blank gap for frame
ex.blbg.ih = ex.rawGaborHeight*2+ex.stim.gapSize;
ex.blbg.iw = ex.rawProbeWidth;

ex.test.blankbg = ones(ex.stim.gapSize+2*ex.rawGaborHeight,ex.rawGaborWidth).*ex.stim.backgroundLum(1,1);
ex.test.blankbgID = Screen('MakeTexture', w, ex.test.blankbg);

%% Sine wave gratings locations 
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+e.vertOffset;
%% bg
xbg = xc;%-ex.stim.horiDistFromFix;
ybg = yc;

%left grating
xL = xc; 
% grating y locations
yL = yc-ex.stim.distFromFix;
%right grating
xR = xc; 
yR = yc+ex.stim.distFromFix; 

%% Fixation cross coordinates
%%%horizontal line
xhorilinel = xc - ex.horiLineW/2;
xhoriliner = xc + ex.horiLineW/2;
yhoriline = yc;

%%%vertical line
xline = xc;
ylineBot = yc + ex.lineH/2;
ylineTop = yc - ex.lineH/2;


%% %%%% initial window - wait for backtick
DrawFormattedText(w,'Fixate the fixation dot as best as you can. \n\n After each drifting stimulus disappears, \n\n 1: Report the MAE effect direction as soon as you see it on the test stimulus \n\n  by pressing the left arrow if central test stimulus appeared to move leftward, \n\n the right arrow if central test stimulus appeared to move rightward, \n\n or down arrow if there was no central MAE \n\n 2: Release the key as soon as central MAE ends \n\n Press Space to start'... % :  '...
    ,xc/2, yc/2,[0 0 0]);
Screen(w, 'Flip', 0);
%WaitSecs(2);
KbTriggerWait(KbName('Space'), deviceNumber);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% START task TASK/FLIPPING
clear l n
n = 1;
blockCnt = 1;
cnt = 0; %
secs = 0;

ex.responses = [];
ex.pressTimes = [];
ex.releaseTimes = [];
ex.resp = [];

% onOffs = [diff([0 ex.longFormBlocks])];
% bLength = ex.blockLength(1);
ex.flipTime = nan(length(ex.trialFlips),length(ex.condShuffle));

KbQueueCreate(deviceNumber,responseKeys);
%%% initial fixation
if n == 1 && blockCnt == 1 %for first block
    ex.tasktstart = clock;
    ex.startRun = GetSecs();
    Screen('FillRect', w, ex.backgroundColor);
    Screen('DrawLines', w, [xline, xline; ylineBot, ylineTop], ex.lineW, [30 30 30]);
    Screen('DrawLines', w, [xhorilinel, xhoriliner; yhoriline, yhoriline],ex.lineW, [30 30 30]);
    Screen(w, 'Flip', 0);
    WaitSecs(ex.initialFixation);
end
%%% Launch the task
for c = 1:length(ex.condShuffle)
    
    condNum = ex.condShuffle(c);
    condName = conditions(condNum).name{:};
    f = 1;
    l = ceil(condNum/2); % indices for each condition
    %l = condNum;
    %flip through the block and following between block time
    while n <= length(ex.longFormBlocks) %true
        KbQueueStart(); % response time
        ex.longFormBlocks(n)
        %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Screen('FillRect', w, ex.stim.backgroundLum);
        if nnz(find(ex.longFormStimOnSecs(n)))
            
            %Background Coordinates
            ex.bgRect = CenterRectOnPoint([0 0 ex.bg.iw ex.bg.ih], xc, yc);
            % Draw bg
            Screen('DrawTexture', w, ex.backgroundID,[], ex.bgRect);
            % Blank gap
            ex.rectCRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.blbg.ih],xc, yc);
            Screen('DrawTexture', w, ex.blankGapID,[],ex.rectCRect);
            ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xL,yL);
            ex.rectRRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xR,yR);

            if contains(condName, 'Sim') && (fix((n-mod(n,ex.flipsPerSec))/(ex.flipsPerSec*2)) - (n-mod(n,ex.flipsPerSec))/(ex.flipsPerSec*2) == 0)


                % stim
                if contains(condName, 'Right')
                    Screen('DrawTexture', w, ex.rectSWaveID(f,l),[],ex.rectLRect);
                    Screen('DrawTexture', w, ex.rectSWaveID(f,l),[],ex.rectRRect);
                    
                elseif contains(condName, 'Left')
                    Screen('DrawTexture', w, ex.rectSWaveID(end-(f-1),l),[],ex.rectLRect);
                    Screen('DrawTexture', w, ex.rectSWaveID(end-(f-1),l),[],ex.rectRRect);
                    
                end
            elseif contains(condName, 'Alt')

                 
                % stim
                if contains(condName, 'Right')
                    if (fix((n-mod(n,ex.flipsPerSec))/(ex.flipsPerSec*2)) - (n-mod(n,ex.flipsPerSec))/(ex.flipsPerSec*2) == 0)
                        Screen('DrawTexture', w, ex.rectSWaveID(f,l),[],ex.rectLRect);
                    else
                        Screen('DrawTexture', w, ex.rectSWaveID(f,l),[],ex.rectRRect);
                    end

                elseif contains(condName, 'Left')
                    if (fix((n-mod(n,ex.flipsPerSec))/(ex.flipsPerSec*2)) - (n-mod(n,ex.flipsPerSec))/(ex.flipsPerSec*2) == 0)

                        Screen('DrawTexture', w, ex.rectSWaveID(end-(f-1),l),[],ex.rectLRect);
                    else
                        Screen('DrawTexture', w, ex.rectSWaveID(end-(f-1),l),[],ex.rectRRect);
                    end
                end

            end
            
        end
        %         Screen('DrawDots', w, [xc yc], ex.fixSize, [255 255 255], [], 2);
        Screen('DrawLines', w, [xline, xline; ylineBot, ylineTop], ex.lineW, [30 30 30]);
        Screen('DrawLines', w, [xhorilinel, xhoriliner; yhoriline, yhoriline],ex.lineW, [30 30 30]);
        
        %% Draw Test stimulus on the screen
        if length(ex.longFormBlocks(1:n))/60 >= ex.blockLength && cnt == 0 %&& length(ex.longFormBlocks(1:n))/60 < ex.blockLength+ex.testLength%(cnt/2 == 1 && GetSecs-time >= 0) && c ~= length(ex.condShuffle)
            %Background Coordinates
            ex.bgRect = CenterRectOnPoint([0 0 ex.bg.iw ex.bg.ih], xc, yc);
            % Draw bg
            Screen('DrawTexture', w, ex.backgroundID,[], ex.bgRect);
            ex.blBgRect = CenterRectOnPoint([0 0 ex.blbg.iw ex.blbg.ih], xc, yc);
            % Draw blank bg
            Screen('DrawTexture', w, ex.test.blankbgID,[], ex.blBgRect);
            
            ex.lcLRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc,yc);
            % test stim
            Screen('DrawTexture', w, ex.testStimID(1),[],ex.lcLRect,[],[],[]);
            %Fixation
            %Screen('DrawDots', w, [xc yc+ex.yoffset], ex.fixSize, [255 255 255], [], 2);
            Screen('DrawLines', w, [xline, xline; ylineBot, ylineTop], ex.lineW, [30 30 30]);
            Screen('DrawLines', w, [xhorilinel, xhoriliner; yhoriline, yhoriline],ex.lineW, [30 30 30]);
            
        end
        
        %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n == 1
            [VBLT, ex.startTrial, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
            flipTimes = ex.startTrial;
            
        else
            [VBLT,flipTime, FlipT, missed] = Screen(w, 'Flip',ex.startTrial + ex.trialFlips(n) - slack); %,   %%% ex.flipTime(n,c)
            flipTimes = [flipTimes, flipTime];
            
        end
        if length(ex.longFormBlocks(1:n))/60>= ex.blockLength
            [secs,~,~] =KbWait(deviceNumber,3);
        end
        if secs %length(ex.longFormBlocks(1:n))/60== ex.blockLength+ex.testLength %&& 
            
            cnt = cnt+1;

            if c ~= length(ex.condShuffle)
                if mod(c,10)
                    DrawFormattedText(w,'Press Space whenever \n\n you feel ready \n\n for the next trial',(4/5)*xc, yc/2,[0 0 0]); %left instruction
                else
                    DrawFormattedText(w,sprintf('Press Space whenever \n\n you feel ready \n\n for trial %d / %d', c+1, length(ex.condShuffle)),(4/5)*xc, yc/2,[0 0 0]); %left instruction
                end
                %% Fixation
                %                 Screen('DrawDots', w, [xc yc], ex.fixSize, [255 255 255], [], 2);
                Screen('DrawLines', w, [xline, xline; ylineBot, ylineTop], ex.lineW, [30 30 30]);
                Screen('DrawLines', w, [xhorilinel, xhoriliner; yhoriline, yhoriline],ex.lineW, [30 30 30]);
                
                Screen(w, 'Flip', 0);
                [~,~,~] = KbWait(deviceNumber,2);
                Screen('FillRect', w, ex.stim.backgroundLum);
                %                 Screen('DrawDots', w, [xc yc], ex.fixSize, [255 255 255], [], 2);
                Screen('DrawLines', w, [xline, xline; ylineBot, ylineTop], ex.lineW, [30 30 30]);
                Screen('DrawLines', w, [xhorilinel, xhoriliner; yhoriline, yhoriline],ex.lineW, [30 30 30]);
                
                Screen(w, 'Flip', 0);
                WaitSecs(ex.ITI2);
                %
            end
        end
        %         if mod(cnt,10) == 0 && cnt >=1
        
        if cnt >=1
            cnt = 0;
            secs = 0;
        end
        KbQueueStop();
        [pressed, firstPress, firstRelease]= KbQueueCheck();
        
        if  (pressed == 1) && ((firstPress(KbName('RightArrow')) > 0 || firstPress( KbName('LeftArrow')) > 0)||(firstPress(KbName('DownArrow')) > 0)) %%
            ex.responses = [ex.responses, 1];
            if (firstPress(KbName('RightArrow')) > 0)
                ex.resp = [ex.resp, 1];
                ex.pressTimes = [ex.pressTimes, firstPress(KbName('RightArrow')) - ex.startRun];
                ex.releaseTimes = [ex.releaseTimes, firstRelease(KbName('RightArrow'))- ex.startRun];

            elseif (firstPress(KbName('LeftArrow')) > 0)
                ex.resp = [ex.resp, 2];
                ex.pressTimes = [ex.pressTimes, firstPress(KbName('LeftArrow')) - ex.startRun];
                ex.releaseTimes = [ex.releaseTimes, firstRelease(KbName('LeftArrow'))- ex.startRun ];

            elseif (firstPress(KbName('DownArrow')) > 0)
                ex.resp = [ex.resp, 3];
                ex.pressTimes = [ex.pressTimes, firstPress(KbName('DownArrow')) - ex.startRun];
                ex.releaseTimes = [ex.releaseTimes, firstRelease(KbName('DownArrow')) - ex.startRun];

            end
               pressed = 0;
        end

        %%%% refresh queue for next character
        KbQueueFlush();
        f = f+1;
        n = n+1;
        if f == ceil(ex.flipsPerSec/ex.stim.cycPerSec)
            f = 1;
        end
    end
    
    ex.flipTime(:,c) = flipTimes;
    n = 1;
end

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

ex.runTime = GetSecs - ex.startRun;

savedir = fullfile(ex.root,'data/static_MAE',sprintf('%s/%s_%s/',ex.version,subject,ex.version));
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(sprintf('/%s_static_MAE_%s_date%s_fix',subject,ex.version,num2str(ex.date)), '.mat'));
%save(savename,'ex');
save(savename,'ex','-v7.3')

ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;  