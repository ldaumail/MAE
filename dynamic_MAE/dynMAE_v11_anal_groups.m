
%MAE Test stimulus with constant contrast v11

names = {'sub-01', 'sub-02','sub-04','sub-05','sub-06','sub-07', 'sub-08', 'sub-09', 'sub-10'}; %'sub-03'
version = 'v11';
responseType = [1, 2, 3];
respFreq = nan(length(responseType), 12, length(names));
% gapSize = [];
for i =1:length(names)
    name = names{i};
   
    %if contains(name, s_f_names) %check if subject had the slow condition first or the fast one
    folderDir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/data',...
        sprintf('/%s/%s_%s',version,name,version), '/');
    
    fileNames = dir(folderDir);
    fileNames = fileNames(contains({fileNames.name},'MAE'));
    fileDir = fileNames.folder;
    exptDat = load(strcat(fileNames(1).folder,'/',fileNames(1).name));
    exptDat = exptDat.ex;
    conditions = exptDat.conds;
    condNums = exptDat.condShuffle;
    
    
    teststarts = exptDat.flipTime(exptDat.blockLength*exptDat.flipsPerSec,1:end)-exptDat.startRun; %exptDat.blockLength*exptDat.flipsPerSec%trial test start time relative to experiment start time
    testends = exptDat.flipTime(end,1:end)-exptDat.startRun+exptDat.ITI1;  %trial test end time relative to experiment start time
   
    
    respDat = zeros(length(responseType), max(exptDat.repsPerRun), exptDat.numConds); %response time data better organized
    % respFreq = nan(length(responseType), max(exptDat.repsPerRun), exptDat.numConds); %response frequency data better organized
    
    cnt = zeros(length(condNums),1);
    
    for c =1:length(condNums)
        condNum = condNums(c);
        cnt(condNum) = cnt(condNum)+1; %counting rep number
        tstart = teststarts(c);
        tend = testends(c);
        respTimes = exptDat.responseTimes(exptDat.responseTimes > tstart &  exptDat.responseTimes < tend);
        resps = exptDat.correctResp(exptDat.responseTimes > tstart &  exptDat.responseTimes < tend);
        resps
        if ~isempty(resps)
            
            respDat(resps(end), cnt(condNum), condNum) = respDat(resps(end), cnt(condNum), condNum)+1; %resps(end) to take the last response in case the participant corrected a mistake
            
        end
    end
    respFreq(:,:,i) = squeeze(sum(respDat,2));
%     gapSize = [gapSize; exptDat.stim.distFromFixDeg];
end
%% Luminance levels check 
% 
% ex.test.contrastOffset = ex.stim.backgroundLum(:,1)./255;% 
% ex.test.contrast = 0.1;
% ex.test.luminanceRange = 2*ex.test.contrast*ex.test.contrastOffset;%0.1; %linspace(0.01,0.20,10);%[0.05, 0.10, 0.15];                                                 % in %, maybe?? %here the number of stimulus contrast levels is the number of different conditions
% ex.test.contrastMultiplicator = ex.test.luminanceRange/2;  % for sine wave 0.5 = 100% contrast, 0.2 = 40%
% 
% bgLum = exptDat.test.contrastOffset*255;
% minLum = bgLum-exptDat.test.contrastMultiplicator*255;
% maxLum = bgLum+exptDat.test.contrastMultiplicator*255;


%% Look at "same" vs "different" direction percepts

same_diff = zeros(3,length(conditions)/2, length(names));
for i = 1:length(names)
    for c =1:length(conditions)
        for t = 1:size(respFreq,1)
            %same direction
            if contains(conditions(c), 'Up') && t == 1
                condcnt = ceil(c/2);
                same_diff(1,condcnt,i) = same_diff(1,condcnt,i)+ respFreq(t,c,i);
            elseif contains(conditions(c), 'Down') && t == 3
                condcnt = ceil(c/2);
                same_diff(1,condcnt,i) = same_diff(1,condcnt,i)+ respFreq(t,c,i);
                %ambiguous direction
            elseif t == 2
                condcnt = ceil(c/2);
                same_diff(2,condcnt,i) = same_diff(2,condcnt,i)+ respFreq(t,c,i);
                %opposite direction
            elseif contains(conditions(c), 'Up') && t == 3
                condcnt = ceil(c/2);
                same_diff(3,condcnt,i) = same_diff(3,condcnt,i)+ respFreq(t,c,i);
            elseif  contains(conditions(c), 'Down') && t == 1
                condcnt = ceil(c/2);
                same_diff(3,condcnt,i) = same_diff(3,condcnt,i)+ respFreq(t,c,i);
            end
        end
    end
end



yvar = 100*[same_diff(:,1,:)./sum(same_diff(:,1,:)), same_diff(:,3,:)./sum(same_diff(:,3,:)), same_diff(:,5,:)./sum(same_diff(:,5,:)),...
     same_diff(:,2,:)./sum(same_diff(:,2,:)), same_diff(:,4,:)./sum(same_diff(:,4,:)), same_diff(:,6,:)./sum(same_diff(:,6,:)),...
    ];

percepts =  {'Same','Ambiguous','Opposite'};
%   {'LowContPhUp'}    {'LowContPhDown'}    {'LowContPhCtUp'}    {'LowContPhCtDown'}    {'MedContPhUp'}    {'MedContPhDown'}    {'MedContPhCtUp'}    {'MedContPhCtDown'}
%    {'HighContPhUp'}    {'HighContPhDown'}    {'HighContPhCtUp'}    {'HighContPhCtDown'}

avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom'),  ...
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom'), ...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom'),... %
    sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom','Control'), ...
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom','Control'), ...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom','Control'),... %
    };
ylab = {'Proportion of occurence of each percept (%)'};
ylims = [0 110];
barPlot2Group(yvar,avgConditions, percepts, ylab, ylims)
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_same_opposite_%s_%s_sp4_test2point67_group.png', name, version)));

%% plot 2: subjects scatter plots + bar dots grouped by percept
percepts =  {'Same','Ambiguous','Opposite'};

avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom'),  ...
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom'), ...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom'),... %
    sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom','Control'), ...
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom','Control'), ...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom','Control'),... %
    };
ylab = {'Proportion of occurence of each percept (%)'};
yval = permute(yvar, [3 1 2]);
%  groupsDotBarPlotSEM(yval,avgConditions, percepts, ylab)
 groupsBarPlotSEM(yval,avgConditions, percepts, ylab)
 plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_same_opposite_%s_sp4_test2point67_group_subplot.svg', version)));


%% stats
% condNames = {'LowPhantom';'MedPhantom';'HighPhantom';'LowPhantomControl';'MedPhantomControl';'HighPhantomControl'};
% 
% conds = repmat((1:length(condNames))',length(names),1);
% 
% subjectsIdx = [];
% for i =1:length(names)
%     subjectsIdx = [subjectsIdx; repmat(names(i),length(condNames),1)];%repmat((1:length(names))',length(condNames),1);
% end

%% Same

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(1,1,:),yvar(1,4,:));%phantom low vs phantom control low
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(1,2,:),yvar(1,5,:));%phantom med vs phantom control med
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(1,3,:),yvar(1,6,:));%phantom high vs phantom control high
[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(1,1,:),yvar(1,2,:));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(1,2,:),yvar(1,3,:));%phantom med vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(1,1,:),yvar(1,3,:));%phantom low vs phantom high
[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(1,4,:),yvar(1,5,:));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(1,5,:),yvar(1,6,:));%phantom control med vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(1,4,:),yvar(1,6,:));%phantom control low vs phantom control high


%% ambiguous

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(2,1,:),yvar(2,4,:));%phantom low vs phantom control low
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(2,2,:),yvar(2,5,:));%phantom med vs phantom control med
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(2,3,:),yvar(2,6,:));%phantom high vs phantom control high
[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(2,1,:),yvar(2,2,:));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(2,2,:),yvar(2,3,:));%phantom med vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(2,1,:),yvar(2,3,:));%phantom low vs phantom high
[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(2,4,:),yvar(2,5,:));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(2,5,:),yvar(2,6,:));%phantom control med vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(2,4,:),yvar(2,6,:));%phantom control low vs phantom control high

%% opposite


[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(3,1,:),yvar(3,4,:));%phantom low vs phantom control low
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(3,2,:),yvar(3,5,:));%phantom med vs phantom control med
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(3,3,:),yvar(3,6,:));%phantom high vs phantom control high
[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(3,1,:),yvar(3,2,:));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(3,2,:),yvar(3,3,:));%phantom med vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(3,1,:),yvar(3,3,:));%phantom low vs phantom high
[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(3,4,:),yvar(3,5,:));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(3,5,:),yvar(3,6,:));%phantom control med vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(3,4,:),yvar(3,6,:));%phantom control low vs phantom control high


%% Inducer type vs contrast interaction

phcond = {'Phantom','Phantom','Phantom','PhantomControl','PhantomControl','PhantomControl'};
contLevels = {'Low';'Med';'High';'Low';'Med';'High'};

contrasts = [];
phantoms = [];
cont =[1 2 3 1 2 3]';
ph = [1 1 1 2 2 2]';

contrasts = repmat(cont,length(names),1);
phantoms = repmat(ph,length(names),1);

subjectsIdx = [];
for i =1:length(names)
    subjectsIdx = [subjectsIdx, repmat(names(i),1,length(contLevels))];%repmat((1:length(names))',length(condNames),1);
end

%% same
data = reshape(yvar(1,:,:), [size(yvar(1,:,:),2)*size(yvar(1,:,:),3),1]);
tbl = table(subjectsIdx', data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);

%% ambiguous
data = reshape(yvar(2,:,:), [size(yvar(2,:,:),2)*size(yvar(2,:,:),3),1]);
tbl = table(subjectsIdx', data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);

%% opposite
data = reshape(yvar(3,:,:), [size(yvar(3,:,:),2)*size(yvar(3,:,:),3),1]);
tbl = table(subjectsIdx', data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);

%% Percent bias

same_diff = zeros(3,length(conditions)/2, length(names));
for i = 1:length(names)
    for c =1:length(conditions)
        for t = 1:size(respFreq,1)
            %same direction
            if contains(conditions(c), 'Up') && t == 1
                condcnt = ceil(c/2);
                same_diff(1,condcnt,i) = same_diff(1,condcnt,i)+ respFreq(t,c,i);
            elseif contains(conditions(c), 'Down') && t == 3
                condcnt = ceil(c/2);
                same_diff(1,condcnt,i) = same_diff(1,condcnt,i)+ respFreq(t,c,i);
                %ambiguous direction
            elseif t == 2
                condcnt = ceil(c/2);
                same_diff(2,condcnt,i) = same_diff(2,condcnt,i)+ respFreq(t,c,i);
                %opposite direction
            elseif contains(conditions(c), 'Up') && t == 3
                condcnt = ceil(c/2);
                same_diff(3,condcnt,i) = same_diff(3,condcnt,i)+ respFreq(t,c,i);
            elseif  contains(conditions(c), 'Down') && t == 1
                condcnt = ceil(c/2);
                same_diff(3,condcnt,i) = same_diff(3,condcnt,i)+ respFreq(t,c,i);
            end
        end
    end
end

%bias
score=[0; 0.5; 1]; %score each percept type: 0 = same, 0.5 = ambiguous, 1 = opposite

bias = same_diff.*score;

percentBias = nan(size(bias,3),size(bias,2));

for i =1:length(names)
    for c =1:size(bias,2)
        percentBias(i,c) = 100*sum(bias(:,c,i))/sum(same_diff(:,c,i));
    end
end

yvar = percentBias;

avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom'),  ...
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom'), ...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom'),... %
    sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom','Control'), ...
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom','Control'), ...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom','Control'),... %
    };
ylab = {'Percent bias (%)'};
ylims = [50 105];

singleBarPlotSEM3(yvar,avgConditions, {'bias'}, ylab, ylims);

% singleBarDotPlotSEM3(yvar,avgConditions, {'bias'}, ylab, ylims);
 plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('percent_bias_%s_sp5.svg', version)));

%% Inducer type vs contrast interaction

phcond = {'Phantom','Phantom','Phantom','PhantomControl','PhantomControl','PhantomControl'};
contLevels = {'Low';'Med';'High';'Low';'Med';'High'};

contrasts = [];
phantoms = [];
cont =[1 2 3 1 2 3]';
ph = [1 1 1 2 2 2]';
for i =1:length(contLevels)
    contrasts = [contrasts; repmat(cont(i),length(names),1)];
    phantoms = [phantoms; repmat(ph(i),length(names),1)];
end

subjectsIdx = repmat(names,1,length(contLevels));%repmat((1:length(names))',length(condNames),1);


%% Percent bias
data = reshape(yvar(:,:), [size(yvar(:,:),1)*size(yvar(:,:),2),1]);
tbl = table(subjectsIdx', data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);


[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1),yvar(:,4));%phantom low vs phantom control low
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,2),yvar(:,5));%phantom med vs phantom control med
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,3),yvar(:,6));%phantom high vs phantom control high
[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,1),yvar(:,2));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,2),yvar(:,3));%phantom med vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,1),yvar(:,3));%phantom low vs phantom high
[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,4),yvar(:,5));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,5),yvar(:,6));%phantom control med vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,4),yvar(:,6));%phantom control low vs phantom control high



%% Correlation between perceptual report and percent bias

%x axis = No | Yes
%Y axis = Percent bias


%No = 0 %Yes =1 % 0.5 = partially, does not extend all the way, but >= 15%

percept = [1 1 0; 0 1 0; 0 1 0; 0.5 1 0; 1 1 1; 1 1 0; 1 1 0; 0 0.5 0; 0 0 0.5]*100;

xval = percept;
yval = yvar(:,1:3);

avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom'),  ...
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom'), ...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom'),... %
       };

xlab = {'Percept strength (%)'};
ylab = {'Percent bias (%)'};
ylims = [0 110];
titl = 'Percent bias as a function of phantom percept';
linregdotsubPlot2(xval,yval, avgConditions, xlab, ylab, titl, ylims)
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, 'percent_bias_phantom_strength_sp5.svg'));

r1 = corrcoef(xval(:,1),yval(:,1));
r2 = corrcoef(xval(:,2),yval(:,2));
r3 = corrcoef(xval(:,3),yval(:,3));

