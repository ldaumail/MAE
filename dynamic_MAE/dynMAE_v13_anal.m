%MAE Test stimulus with constant background luminance v13

names = {'sub-01','sub-02','sub-03', 'sub-04', 'sub-05', 'sub-06', 'sub-07','sub-08','sub-09', 'sub-10', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15'};%
version = 'v13';
responseType = [1, 2, 3];
respFreq = nan(length(responseType), 18, length(names));

for i =1:length(names)
    name = names{i};
   
    %if contains(name, s_f_names) %check if subject had the slow condition first or the fast one
    folderDir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/data',...
        sprintf('/dynamic_MAE/%s/%s_%s',version,name,version), '/');
    
    fileNames = dir(folderDir);
    fileNames = fileNames(contains({fileNames.name},'MAE'));
    fileDir = fileNames.folder;
    exptDat = load(strcat(fileNames(1).folder,'/',fileNames(1).name));
    exptDat = exptDat.ex;
    conditions1 = exptDat.conds;
    condNums = exptDat.condShuffle;
    
    
    teststarts = exptDat.flipTime(exptDat.blockLength*exptDat.flipsPerSec,1:end)-exptDat.startRun; %exptDat.blockLength*exptDat.flipsPerSec%trial test start time relative to experiment start time
    testends = exptDat.flipTime(end,1:end)-exptDat.startRun+exptDat.ITI1;  %trial test end time relative to experiment start time
   
    
    respDat = zeros(length(responseType), max(exptDat.repsPerRun), exptDat.numConds); %response time data better organized
    
    cnt = zeros(length(condNums),1);
    
    for c =1:length(condNums)
        condNum = condNums(c);
        cnt(condNum) = cnt(condNum)+1; %counting rep number
        tstart = teststarts(c);
        tend = testends(c);
        respTimes = exptDat.responseTimes(exptDat.responseTimes > tstart &  exptDat.responseTimes < tend);
        resps = exptDat.resp(exptDat.responseTimes > tstart &  exptDat.responseTimes < tend);
        resps
        if ~isempty(resps)
            
            respDat(resps(end), cnt(condNum), condNum) = respDat(resps(end), cnt(condNum), condNum)+1; %resps(end) to take the last response in case the participant corrected a mistake
            
        end
    end
    respFreq(:,1:length(conditions1),i) = squeeze(sum(respDat,2));

    clear exptDat condNums teststarts testends cnt
    exptDat = load(strcat(fileNames(2).folder,'/',fileNames(2).name));
    exptDat = exptDat.ex;
    conditions2 = exptDat.conds;
    condNums = exptDat.condShuffle;
    
    
    teststarts = exptDat.flipTime(exptDat.blockLength*exptDat.flipsPerSec,1:end)-exptDat.startRun; %exptDat.blockLength*exptDat.flipsPerSec%trial test start time relative to experiment start time
    testends = exptDat.flipTime(end,1:end)-exptDat.startRun+exptDat.ITI1;  %trial test end time relative to experiment start time
   
    
    respDat = zeros(length(responseType), max(exptDat.repsPerRun), exptDat.numConds); %response time data better organized
    
    cnt = zeros(length(condNums),1);
    
    for c =1:length(condNums)
        condNum = condNums(c);
        cnt(condNum) = cnt(condNum)+1; %counting rep number
        tstart = teststarts(c);
        tend = testends(c);
        respTimes = exptDat.responseTimes(exptDat.responseTimes > tstart &  exptDat.responseTimes < tend);
        resps = exptDat.resp(exptDat.responseTimes > tstart &  exptDat.responseTimes < tend);
        resps
        if ~isempty(resps)
            
            respDat(resps(end), cnt(condNum), condNum) = respDat(resps(end), cnt(condNum), condNum)+1; %resps(end) to take the last response in case the participant corrected a mistake
            
        end
    end
    respFreq(:,length(conditions1)+1:length(conditions1)+length(conditions2),i) = squeeze(sum(respDat,2));

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
conditions = {conditions1{:}, conditions2{:}};
same_diff = zeros(length(names),length(conditions)/2, 3);
for i = 1:length(names)
    for c =1:length(conditions)
        for t = 1:size(respFreq,1)
            %same direction
            if contains(conditions(c), 'Right') && t == 1
                condcnt = ceil(c/2);
                same_diff(i,condcnt,1) = same_diff(i,condcnt,1)+ respFreq(t,c,i);
            elseif contains(conditions(c), 'Left') && t == 2
                condcnt = ceil(c/2);
                same_diff(i,condcnt,1) = same_diff(i,condcnt,1)+ respFreq(t,c,i);
                %ambiguous direction
            elseif t == 3
                condcnt = ceil(c/2);
                same_diff(i,condcnt,2) = same_diff(i,condcnt,2)+ respFreq(t,c,i);
                %opposite direction
            elseif contains(conditions(c), 'Right') && t == 2
                condcnt = ceil(c/2);
                same_diff(i,condcnt,3) = same_diff(i,condcnt,3)+ respFreq(t,c,i);
            elseif  contains(conditions(c), 'Left') && t == 1
                condcnt = ceil(c/2);
                same_diff(i,condcnt,3) = same_diff(i,condcnt,3)+ respFreq(t,c,i);
            end
        end
    end
end

yvar= [];
for c = 1:size(same_diff,2)
    yvar = [yvar, 100*same_diff(:,c,:)./sum(same_diff(:,c,:),3)];
end
% percepts =  {'Same','Ambiguous','Opposite'};
% 
% avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Low Contrast Full'),  ...
%     sprintf('%s\\newline%s\\newline%s\n','Med Contrast Full'), ...
%     sprintf('%s\\newline%s\\newline%s\n','High Contrast Full'),... %
%     sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom'),  ...
%     sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom Control'), ...
%     sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom'),... 
%     sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom Control'),... 
%     sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom'),... 
%     sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom Control'),... 
%     };
% ylab = {'Proportion of occurence of each percept (%)'};
% ylims = [0 110];
% % barPlot2Group(yvar,avgConditions, percepts, ylab, ylims)
% singleBarPlot(yvar(:,:,2),avgConditions, percepts, ylab, ylims)
% % singleBarPlotSEM2(yvar,avgConditions, percepts, ylab, ylims)
% plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
% mkdir(plotdir);
% saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_%s_%s.png', name, version)));

%% plot 2: subjects scatter plots + bar dots grouped by percept
percepts =  {'Same','Ambiguous','Opposite'};

% avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Low Contrast Full'),  ...
%     sprintf('%s\\newline%s\\newline%s\n','Med Contrast Full'), ...
%     sprintf('%s\\newline%s\\newline%s\n','High Contrast Full'),... %
%     sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom'),  ...
%     sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom Control'), ...
%     sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom'),... 
%     sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom Control'),... 
%     sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom'),... 
%     sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom Control'),... 
%     };
ylab = {'Proportion of occurence of each percept (%)'};
% yval = permute(yvar, [1 3 2]); %respType * subjs * all conds
%%groupsDotBarPlotSEM(yval,avgConditions, percepts, ylab)
%   groupsBarPlotSEM(yval,avgConditions, percepts, ylab)

 
contLevel = {'Low', 'Medium', 'High'};
conds = {'Full Grating','Phantom','Phantom Control'};
% condMAEDir is subjs * all conds * respType here
condMAEDir = [yvar(:,1,:), yvar(:,2,:), yvar(:,3,:), yvar(:,4,:), yvar(:,6,:), yvar(:,8,:),  yvar(:,5,:), yvar(:,7,:), yvar(:,9,:)];

yvar2 = reshape(condMAEDir, length(names), length(contLevel), length(conds), length(responseType));
ylims = [0 20; 0 20; 60 105];
groupsBarPlotSEM2(yvar2,conds, contLevel, percepts, ylab,ylims)
 plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
 saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_same_opposite_%s_group_subplot.svg', version)));


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

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,4,1),yvar(:,5,1));%phantom low vs phantom control low
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,6,1),yvar(:,7,1));%phantom med vs phantom control med
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,8,1),yvar(:,9,1));%phantom high vs phantom control high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,1,1),yvar(:,4,1));%full low vs phantom low 
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,2,1),yvar(:,6,1));%full med vs phantom med 
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,3,1),yvar(:,8,1));%full high vs phantom high

[ttestMean(4), Pval(7),~,Stats(7).stats] = ttest(yvar(:,1,1),yvar(:,5,1));%full low vs phantom control low
[ttestMean(5), Pval(8),~,Stats(8).stats] = ttest(yvar(:,2,1),yvar(:,7,1));%full med vs phantom control med
[ttestMean(6), Pval(9),~,Stats(9).stats] = ttest(yvar(:,3,1),yvar(:,9,1));%full high vs phantom control high




%% ambiguous

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,4,2),yvar(:,5,2));%phantom low vs phantom control low
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,6,2),yvar(:,7,2));%phantom med vs phantom control med
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,8,2),yvar(:,9,2));%phantom high vs phantom control high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,1,2),yvar(:,4,2));%full low vs phantom low 
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,2,2),yvar(:,6,2));%full med vs phantom med 
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,3,2),yvar(:,8,2));%full high vs phantom high

[ttestMean(4), Pval(7),~,Stats(7).stats] = ttest(yvar(:,1,2),yvar(:,5,2));%full low vs phantom control low
[ttestMean(5), Pval(8),~,Stats(8).stats] = ttest(yvar(:,2,2),yvar(:,7,2));%full med vs phantom control med
[ttestMean(6), Pval(9),~,Stats(9).stats] = ttest(yvar(:,3,2),yvar(:,9,2));%full high vs phantom control high
% 
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
conditions = {conditions1{:}, conditions2{:}};

same_diff = zeros(3,length(conditions)/2, length(names));
for i = 1:length(names)
    for c =1:length(conditions)
        for t = 1:size(respFreq,1)
            %same direction
            if contains(conditions(c), 'Right') && t == 1
                condcnt = ceil(c/2);
                same_diff(1,condcnt,i) = same_diff(1,condcnt,i)+ respFreq(t,c,i);
            elseif contains(conditions(c), 'Left') && t == 2
                condcnt = ceil(c/2);
                same_diff(1,condcnt,i) = same_diff(1,condcnt,i)+ respFreq(t,c,i);
                %ambiguous direction
            elseif t == 3
                condcnt = ceil(c/2);
                same_diff(2,condcnt,i) = same_diff(2,condcnt,i)+ respFreq(t,c,i);
                %opposite direction
            elseif contains(conditions(c), 'Right') && t == 2
                condcnt = ceil(c/2);
                same_diff(3,condcnt,i) = same_diff(3,condcnt,i)+ respFreq(t,c,i);
            elseif  contains(conditions(c), 'Left') && t == 1
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


contLevel = {'Low', 'Medium', 'High'};
conds = {'Full Grating','Phantom','Phantom Control'};

condPercentBias = [percentBias(:,1), percentBias(:,2), percentBias(:,3), percentBias(:,4), percentBias(:,6), percentBias(:,8),  percentBias(:,5), percentBias(:,7), percentBias(:,9)];

yvar = reshape(condPercentBias, length(names), length(contLevel),length(conds));

ylab = {'Percent bias (%)'};
ylims = [40 110];

%%singleBarLinePlotSEM(percentBias',avgConditions, ylab, ylims)%  singleBarPlot(yvar(:,1)',avgConditions, {'bias'}, ylab, ylims);
MAEBarPlotSEM(yvar,conds, ylab, ylims, contLevel)

% singleBarDotPlotSEM3(yvar,avgConditions, {'bias'}, ylab, ylims);
 plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('dynamicMAE_percent_bias_%s.png', version)));
%% %% Inducer type vs contrast interaction

phcond = {'Full','Full','Full','Phantom','Phantom','Phantom','PhantomControl','PhantomControl','PhantomControl'};
contLevels = {'Low';'Med';'High';'Low';'Med';'High';'Low';'Med';'High'};

contrasts = [];
phantoms = [];
cont =[1 2 3 1 2 3 1 2 3]';
ph = [1 1 1 2 2 2 3 3 3]';
for i =1:length(contLevels)
    contrasts = [contrasts; repmat(cont(i),length(names),1)];
    phantoms = [phantoms; repmat(ph(i),length(names),1)];
end

subjectsIdx = repmat(names',length(contLevels),1);%repmat((1:length(names))',length(condNames),1);


data = reshape(condPercentBias, [size(condPercentBias,1)*size(condPercentBias,2),1]);
tbl = table(subjectsIdx, data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);

%% ttests
yvar = condPercentBias;

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,4),yvar(:,7));%phantom low vs phantom control low
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,5),yvar(:,8));%phantom med vs phantom control med
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,6),yvar(:,9));%phantom high vs phantom control high
[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4),yvar(:,5));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,5),yvar(:,6));%phantom med vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,4),yvar(:,6));%phantom low vs phantom high
[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,7),yvar(:,8));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,8),yvar(:,9));%phantom control med vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,6),yvar(:,9));%phantom control low vs phantom control high


%% Histogram plot

yvar = percentBias;
avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Low Contrast Full'),  ...
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Full'), ...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast Full'),... %
    sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom'),  ...
    sprintf('%s\\newline%s\\newline%s\n','Low Contrast Phantom','Control'), ...
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom'),... %
    sprintf('%s\\newline%s\\newline%s\n','Med Contrast Phantom','Control'), ...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom'), ...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast Phantom','Control'),... %
    };
conds = {'Phantom','Phantom control'};
xlab = {'Percent bias (%)'};
ylab = {'Frequency'};
percentHistogram(yvar',avgConditions, xlab, ylab,[0 5])
 plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('percent_bias_histogram_%s_sp5.png', version)));

% %% Inducer type vs contrast interaction
% 
% phcond = {'Phantom','Phantom','Phantom','PhantomControl','PhantomControl','PhantomControl'};
% contLevels = {'Low';'Med';'High';'Low';'Med';'High'};
% 
% contrasts = [];
% phantoms = [];
% cont =[1 2 3 1 2 3]';
% ph = [1 1 1 2 2 2]';
% for i =1:length(contLevels)
%     contrasts = [contrasts; repmat(cont(i),length(names),1)];
%     phantoms = [phantoms; repmat(ph(i),length(names),1)];
% end
% 
% subjectsIdx = repmat(names,1,length(contLevels));%repmat((1:length(names))',length(condNames),1);
% 

%% Percent bias
% data = reshape(yvar(:,:), [size(yvar(:,:),1)*size(yvar(:,:),2),1]);
% tbl = table(subjectsIdx', data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
% lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %
% 
% [pVal, F, R] = coefTest(lme);


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

