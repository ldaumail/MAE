
names = {'sub-01', 'sub-02','sub-03','sub-04','sub-05','sub-06','sub-08'}; %
version = 'v9';
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
 
%% Look at frequency of the responses for each condition

% yvar = 100*[respFreq(:,1)./sum(respFreq(:,1)), respFreq(:,2)./sum(respFreq(:,2)), respFreq(:,5)./sum(respFreq(:,5)), respFreq(:,6)./sum(respFreq(:,6)), respFreq(:,3)./sum(respFreq(:,3)), respFreq(:,4)./sum(respFreq(:,4))];
yvar = 100*[respFreq(:,1)./sum(respFreq(:,1)), respFreq(:,2)./sum(respFreq(:,2)), respFreq(:,3)./sum(respFreq(:,3)),...
    respFreq(:,4)./sum(respFreq(:,4)), respFreq(:,5)./sum(respFreq(:,5)), respFreq(:,6)./sum(respFreq(:,6)),...
    respFreq(:,7)./sum(respFreq(:,7)), respFreq(:,8)./sum(respFreq(:,8)), respFreq(:,9)./sum(respFreq(:,9)),...
    respFreq(:,10)./sum(respFreq(:,10)), respFreq(:,11)./sum(respFreq(:,11)), respFreq(:,12)./sum(respFreq(:,12)),...
    ];

percepts =  {'Probe going up','Counterphasing','Probe going down'};
%   {'LowContPhUp'}    {'LowContPhDown'}    {'LowContPhCtUp'}    {'LowContPhCtDown'}    {'MedContPhUp'}    {'MedContPhDown'}    {'MedContPhCtUp'}    {'MedContPhCtDown'}
%    {'HighContPhUp'}    {'HighContPhDown'}    {'HighContPhCtUp'}    {'HighContPhCtDown'}

avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom',' up'), sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom',' down'), ...
    sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom','Control up'), sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom','Control down'), ...
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom',' up'),sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom',' down'), ...
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom','Control up'),sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom','Control down'), ...
   sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom',' up'),sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom',' down'),... % 
      sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom','Control up'),sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom','Control down'),... %
 }; 
ylab = {'Proportion of occurence of each percept (%)'};
ylims = [0 110];
barPlot(yvar,avgConditions, percepts, ylab, ylims)
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
% saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_across_conditions_%s_%s_sp4_test2point67.png', name, version)));


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
groupsDotBarPlotSEM(yval,avgConditions, percepts, ylab)
saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_same_opposite_%s_%s_sp4_test2point67_group_subplot.png', name, version)));


%% stats
condNames = {'LowPhantom';'MedPhantom';'HighPhantom';'LowPhantomControl';'MedPhantomControl';'HighPhantomControl'};

conds = repmat((1:length(condNames))',length(names),1);

subjectsIdx = [];
for i =1:length(names)
    subjectsIdx = [subjectsIdx; repmat(names(i),length(condNames),1)];%repmat((1:length(names))',length(condNames),1);
end

%% Same
data = reshape(yvar(1,:,:), [size(yvar(1,:,:),2)*size(yvar(1,:,:),3),1]);
tbl = table(subjectsIdx, data, conds,'VariableNames',{'SubjectIndex','ResponseFreq','Condition'});
lme = fitlme(tbl,'ResponseFreq~Condition+(1|SubjectIndex)+(Condition-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);

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
data = reshape(yvar(2,:,:), [size(yvar(2,:,:),2)*size(yvar(2,:,:),3),1]);
tbl = table(subjectsIdx, data, conds,'VariableNames',{'SubjectIndex','ResponseFreq','Condition'});
lme = fitlme(tbl,'ResponseFreq~Condition+(1|SubjectIndex)+(Condition-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);

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

data = reshape(yvar(3,:,:), [size(yvar(3,:,:),2)*size(yvar(3,:,:),3),1]);
tbl = table(subjectsIdx, data, conds,'VariableNames',{'SubjectIndex','ResponseFreq','Condition'});
lme = fitlme(tbl,'ResponseFreq~Condition+(1|SubjectIndex)+(Condition-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);

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
ylims = [50 110];
% groupsBarPlotSEM(yvar',avgConditions, {'bias'}, ylab)
% groupsDotBarPlotSEM(yvar',avgConditions, {'bias'}, ylab)

% singleBarPlotSEM3(yvar,avgConditions, {'bias'}, ylab, ylims);

singleBarDotPlotSEM3(yvar,avgConditions, {'bias'}, ylab, ylims);
saveas(gcf,strcat(plotdir, sprintf('percent_bias_%s_%s_sp4_test2point67.png', name, version)));

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


