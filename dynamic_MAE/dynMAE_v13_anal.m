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


%% Look at "same" vs "ambiguous" vs "different" direction percepts
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

%% plot 2: subjects scatter plots + bar dots grouped by percept
percepts =  flip({'Same','Ambiguous','Opposite'});

ylab = {'Proportion of occurence of each percept (%)'};

contLevel = {'Low', 'Medium', 'High'};
conds = {'Full Grating','Phantom','Phantom Control'};
% condMAEDir is subjs * all conds * respType here
condMAEDir = [yvar(:,1,:), yvar(:,2,:), yvar(:,3,:), yvar(:,4,:), yvar(:,6,:), yvar(:,8,:),  yvar(:,5,:), yvar(:,7,:), yvar(:,9,:)];

yvar2 = flip(reshape(condMAEDir, length(names), length(contLevel), length(conds), length(responseType)),4);
ylims = flip([0 50; 0 50; 50 100]);
groupsBarPlotSEM2(yvar2,conds, contLevel, percepts, ylab,ylims)
 plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
 saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_same_opposite_%s_group_subplot.svg', version)));


%% stats

%% Same

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1,1),yvar(:,2,1));%full low vs full med 
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,1,1),yvar(:,3,1));%full low vs full high 
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,2,1),yvar(:,3,1));%full med vs full high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4,1),yvar(:,6,1));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,4,1),yvar(:,8,1));%phantom low vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,6,1),yvar(:,8,1));%phantom med vs phantom high

[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,5,1),yvar(:,7,1));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,5,1),yvar(:,9,1));%phantom control low vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,7,1),yvar(:,9,1));%phantom control med vs phantom control high




%% ambiguous

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1,2),yvar(:,2,2));%full low vs full med 
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,1,2),yvar(:,3,2));%full low vs full high 
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,2,2),yvar(:,3,2));%full med vs full high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4,2),yvar(:,6,2));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,4,2),yvar(:,8,2));%phantom low vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,6,2),yvar(:,8,2));%phantom med vs phantom high

[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,5,2),yvar(:,7,2));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,5,2),yvar(:,9,2));%phantom control low vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,7,2),yvar(:,9,2));%phantom control med vs phantom control high
% 
%% opposite


[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1,3),yvar(:,2,3));%full low vs full med 
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,1,3),yvar(:,3,3));%full low vs full high 
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,2,3),yvar(:,3,3));%full med vs full high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4,3),yvar(:,6,3));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,4,3),yvar(:,8,3));%phantom low vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,6,3),yvar(:,8,3));%phantom med vs phantom high

[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,5,3),yvar(:,7,3));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,5,3),yvar(:,9,3));%phantom control low vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,7,3),yvar(:,9,3));%phantom control med vs phantom control high


%% Inducer type vs contrast interaction

phcond = {'Full','Full','Full','Phantom','Phantom','Phantom','PhantomControl','PhantomControl','PhantomControl'};
contLevels = {'Low';'Med';'High';'Low';'Med';'High';'Low';'Med';'High'};


ph = [1 1 1 2 2 2 3 3 3]';
cont =[1 2 3 1 2 3 1 2 3]';

contrasts = [];
phantoms = [];
for c =1:length(cont)
    contrasts = [contrasts; repmat(cont(c),length(names),1)];
    phantoms = [phantoms; repmat(ph(c),length(names),1)];
end
subjectsIdx = [];
subjectsIdx = repmat(names,1,length(contLevels));%repmat((1:length(names))',length(condNames),1);

% %% same
data = reshape(condMAEDir(:,:,1), [size(condMAEDir(:,:,1),1)*size(condMAEDir(:,:,1),2),1]);
tbl = table(subjectsIdx', data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);

% 
%% ambiguous
data = reshape(condMAEDir(:,:,2), [size(condMAEDir(:,:,2),1)*size(condMAEDir(:,:,2),2),1]);
tbl = table(subjectsIdx', data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);
% 
% %% opposite
data = reshape(condMAEDir(:,:,3), [size(condMAEDir(:,:,3),1)*size(condMAEDir(:,:,3),2),1]);
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

%% Plot
contLevel = {'Low', 'Medium', 'High'};
conds = {'Full Grating','Phantom','Phantom Control'};

condPercentBias = [percentBias(:,1), percentBias(:,2), percentBias(:,3), percentBias(:,4), percentBias(:,6), percentBias(:,8),  percentBias(:,5), percentBias(:,7), percentBias(:,9)];

yvar = reshape(condPercentBias, length(names), length(contLevel),length(conds));

ylab = {'Percent bias (%)'};
ylims = [50 110];

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
ph = [1 1 1 2 2 2 3 3 3]';
cont =[1 2 3 1 2 3 1 2 3]';



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

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1),yvar(:,2));%full low vs full med 
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,1),yvar(:,3));%full low vs full high 
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,2),yvar(:,3));%full med vs full high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4),yvar(:,5));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,4),yvar(:,6));%phantom low vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,5),yvar(:,6));%phantom med vs phantom high

[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,7),yvar(:,8));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,7),yvar(:,9));%phantom control low vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,8),yvar(:,9));%phantom control med vs phantom control high


%% Subtract Same from Opposite Percent Bias

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
score=[-1; 0; 1]; %score each percept type: 0 = same, 0.5 = ambiguous, 1 = opposite

bias = same_diff.*score;

percentBias = nan(size(bias,3),size(bias,2));

for i =1:length(names)
    for c =1:size(bias,2)
        percentBias(i,c) = 100*sum(bias(:,c,i))/sum(same_diff(:,c,i));
    end
end

%% Plot
contLevel = {'Low', 'Medium', 'High'};
conds = {'Full Grating','Phantom','Phantom Control'};

condPercentBias = [percentBias(:,1), percentBias(:,2), percentBias(:,3), percentBias(:,4), percentBias(:,6), percentBias(:,8),  percentBias(:,5), percentBias(:,7), percentBias(:,9)];

yvar = reshape(condPercentBias, length(names), length(contLevel),length(conds));

ylab = {'Percent bias (%)'};
ylims = [50 100];

%%singleBarLinePlotSEM(percentBias',avgConditions, ylab, ylims)%  singleBarPlot(yvar(:,1)',avgConditions, {'bias'}, ylab, ylims);
MAEBarPlotSEM(yvar,conds, ylab, ylims, contLevel)

% singleBarDotPlotSEM3(yvar,avgConditions, {'bias'}, ylab, ylims);
 plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('dynamicMAE_subtract_percent_bias_%s.svg', version)));

%% stats

yvar = condPercentBias;

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1),yvar(:,2));%full low vs full med 
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,1),yvar(:,3));%full low vs full high 
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,2),yvar(:,3));%full med vs full high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4),yvar(:,5));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,4),yvar(:,6));%phantom low vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,5),yvar(:,6));%phantom med vs phantom high

[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,7),yvar(:,8));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,7),yvar(:,9));%phantom control low vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,8),yvar(:,9));%phantom control med vs phantom control high

%% Inducer type vs contrast interaction

phcond = {'Full','Full','Full','Phantom','Phantom','Phantom','PhantomControl','PhantomControl','PhantomControl'};
contLevels = {'Low';'Med';'High';'Low';'Med';'High';'Low';'Med';'High'};

contrasts = [];
phantoms = [];
ph = [1 1 1 2 2 2 3 3 3]';
cont =[1 2 3 1 2 3 1 2 3]';



for i =1:length(contLevels)
    contrasts = [contrasts; repmat(cont(i),length(names),1)];
    phantoms = [phantoms; repmat(ph(i),length(names),1)];
end

subjectsIdx = repmat(names',length(contLevels),1);%repmat((1:length(names))',length(condNames),1);


data = reshape(condPercentBias, [size(condPercentBias,1)*size(condPercentBias,2),1]);
tbl = table(subjectsIdx, data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);

%%
% 
% phResps = nan(2,12,length(names));
% for i =1:length(names)
%     name = names{i};
%     folderDir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/data',...
%         sprintf('/dynamic_MAE/%s/%s_%s',version,name,version), '/');
%     fileNames = dir(folderDir);
%     fileNames = string({fileNames(contains({fileNames.name},'percept')).name});
% 
%     exptName = strcat(folderDir,fileNames(contains(fileNames,'percept')));
%     exptDat = load(exptName);
%     exptDat = exptDat.ex;
%     conditions = exptDat.conds;
%     condNums = exptDat.condShuffle;
%     
%     trConds = exptDat.condResp; %time at which a cycle of a condition starts. May be multiple cycles per trial for allowing time for participant to respond
%     trChange = [1, diff(trConds)]; %Focus on cycles for which the condition changed = trial starts
%     cycleStarts = exptDat.trialStartTime-exptDat.startRun; %Time that a cycle starts %exptDat.flipTime(exptDat.blockLength*exptDat.flipsPerSec,1:end)-exptDat.startRun; %exptDat.blockLength*exptDat.flipsPerSec%trial test start time relative to experiment start time
%     trstarts = cycleStarts(cycleStarts & trChange); %Get trial start time
%     trRespsT = exptDat.responseTimes; % Get response time
%     trResps = exptDat.resp; %Get response type
%     
%     cnt = zeros(length(unique(condNums)),1);
%     c = 0;
%     for t =1:length(trstarts)
%         condNum = condNums(t);
%         cnt(condNum) = cnt(condNum)+1; %counting rep number
%         if t < length(trstarts)
%             if nnz((trRespsT > trstarts(t) & trRespsT < trstarts(t + 1)))
%                 resp = exptDat.resp((trRespsT > trstarts(t) & trRespsT < trstarts(t + 1)));
%                 resp = resp(end);
%                 phResps(cnt(condNum),condNum,i) = resp;
%             end
%         elseif t == length(trstarts) & nnz((trRespsT > trstarts(t)))
%                 resp = exptDat.resp((trRespsT > trstarts(t)));
%                 resp = resp(end);
%                 phResps(cnt(condNum),condNum,i) = resp;
%         end
%     end
%     
% end
% 
% 
% orgPhResps = squeeze(nanmean([[phResps(:,1,:);phResps(:,2,:)],[phResps(:,5,:);phResps(:,6,:)],[phResps(:,9,:);phResps(:,10,:)],[phResps(:,3,:);phResps(:,4,:)],[phResps(:,7,:);phResps(:,8,:)],[phResps(:,11,:);phResps(:,12,:)]],1))';
% 
% condPhResps = reshape(orgPhResps,size(orgPhResps,1),3, 2); % %subj resp x contrast level x cond 
% 
% condBiasResps = reshape(condPercentBias, length(names), length(contLevel),length(conds));
% 
% 
% xval = reshape(condPhResps, size(condPhResps,1), size(condPhResps,2)*size(condPhResps,3)).*100/5; %
% 
% yval =  reshape(condBiasResps(:,:,2:3), size(condBiasResps(:,:,2:3),1),size(condBiasResps(:,:,2:3),2)*size(condBiasResps(:,:,2:3),3));
% conditions = {'Low Contrast Phantom','Medium Contrast Phantom','High Contrast Phantom', 'Low Contrast Phantom Control','Medium Contrast Phantom Control','High COntrast Phantom Control',}; 
% xlabels = {'Phantom Vividness'};
% ylabels = {'MAE Bias'};
% titl = {'MAE bias as a function of perceived phantom vividness'};
% ylims =[0 100];
% linregdotsubPlot(xval,yval,conditions, xlabels, ylabels, titl, ylims)
% 
%  plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
% mkdir(plotdir);
% saveas(gcf,strcat(plotdir, sprintf('dynamicMAE_percept_strength_vs_percent_bias_%s.svg', version)));
% 
% 
% %% Linear trend of MAE bias as a function of percept strength for a given contrast level
% contLevel = {'Low', 'Medium', 'High'};
% xval = reshape(condPhResps, size(condPhResps,1), size(condPhResps,2)*size(condPhResps,3)).*100/5; %
% yval =  reshape(condBiasResps(:,:,2:3), size(condBiasResps(:,:,2:3),1),size(condBiasResps(:,:,2:3),2)*size(condBiasResps(:,:,2:3),3));
% ylims =[0 100];
% 
% %Phantom
% figure()
% scatter(mean(xval(:,1:length(contLevel)),1),mean(yval(:,1:length(contLevel)),1), 'blue')
% hold on
% plot(mean(xval(:,1:length(contLevel)),1),mean(yval(:,1:length(contLevel)),1),'col', 'blue')
% hold on
% coeffs = polyfit(mean(xval(:,1:length(contLevel)),1),mean(yval(:,1:length(contLevel)),1),1);
% f = polyval(coeffs,mean(xval(:,1:length(contLevel)),1));
% plot(mean(xval(:,1:length(contLevel)),1),mean(yval(:,1:length(contLevel)),1),'o',mean(xval(:,1:length(contLevel)),1), f,'-','Color',[0 0 1],'MarkerSize',2, 'MarkerFaceColor',[0 0 1],'linewidth',2)
%         
% linreg = fitlm(mean(xval(:,1:length(contLevel)),1),mean(yval(:,1:length(contLevel)),1));
% [linPvalue(1,1),F(1,1),r(1,1)] = coefTest(linreg);
% Rsquared(1,1) = linreg.Rsquared.Ordinary; %R2
% pearsonR = corrcoef(mean(xval(:,1:length(contLevel)),1),mean(yval(:,1:length(contLevel)),1));
% t = text(30,25, sprintf('y = %.2f + %.2f*x',round(coeffs(2),2), round(coeffs(1),2)),'fontweight','bold','fontsize',10);
% t(1).Color = [0 0 1];
% t = text(30,20, sprintf('p = %.2f, r^2 = %.2f', linPvalue, Rsquared),'fontweight','bold','fontsize',10);
% t(1).Color = [0 0 1];
% t = text(30,15, sprintf('Pearson r = %.2f', pearsonR(1,2)),'fontweight','bold','fontsize',10);
% t(1).Color = [0 0 1];
% %Phantom control
% scatter(mean(xval(:,length(contLevel)+1:end),1),mean(yval(:,length(contLevel)+1:end),1), 'red')
% hold on
% plot(mean(xval(:,length(contLevel)+1:end),1),mean(yval(:,length(contLevel)+1:end,1)), 'col', 'red')
% hold on
% plot(min(ylims):max(ylims)+1, min(ylims):max(ylims)+1, 'Color', [0.75 0.75 0.75],'linewidth',2)
% 
% coeffs = polyfit(mean(xval(:,length(contLevel)+1:end),1),mean(yval(:,length(contLevel)+1:end,1)),1);
% f = polyval(coeffs,mean(xval(:,length(contLevel)+1:end),1));
% plot(mean(xval(:,length(contLevel)+1:end),1),mean(yval(:,length(contLevel)+1:end,1)),'o',mean(xval(:,length(contLevel)+1:end),1), f,'-','Color',[1 0 0],'MarkerSize',2, 'MarkerFaceColor',[1 0 0],'linewidth',2)
%         
% linreg = fitlm(mean(xval(:,length(contLevel)+1:end),1),mean(yval(:,length(contLevel)+1:end,1)));
% [linPvalue(1,1),F(1,1),r(1,1)] = coefTest(linreg);
% Rsquared(1,1) = linreg.Rsquared.Ordinary; %R2
% pearsonR = corrcoef(mean(xval(:,length(contLevel)+1:end),1),mean(yval(:,length(contLevel)+1:end,1)));
% t = text(30,50, sprintf('y = %.2f + %.2f*x',round(coeffs(2),2), round(coeffs(1),2)),'fontweight','bold','fontsize',10);
% t(1).Color = [1 0 0];
% t = text(30,45, sprintf('p = %.2f, r^2 = %.2f', linPvalue, Rsquared),'fontweight','bold','fontsize',10);
% t(1).Color = [1 0 0];
% t = text(30,40, sprintf('Pearson r = %.2f', pearsonR(1,2)),'fontweight','bold','fontsize',10);
% t(1).Color = [1 0 0];
% xlabel('Phantom Vividness')
% ylabel('MAE bias (%)')
%  plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
% mkdir(plotdir);
% saveas(gcf,strcat(plotdir, sprintf('dynamicMAE_percept_strength_vs_percent_bias_contrast_linearreg_%s.svg', version)));
% 
% %% Make same plot with normalized data
% 
% 
% %% Difference between phantom and phantom control condition
% 
% yval =  reshape(condBiasResps(:,:,2:3), size(condBiasResps(:,:,2:3),1),size(condBiasResps(:,:,2:3),2)*size(condBiasResps(:,:,2:3),3));
% 
% biasDiff = yval(:,1:length(contLevel))-yval(:,length(contLevel)+1:end);
% biasDiff = reshape(biasDiff, [size(biasDiff), 1]);
% ylab ={'Percent bias difference'};
% 
% 
% MAEdiffBarPlotSEM(biasDiff,{''}, ylab, [-10 10], contLevel)
% 
% %% Stats
% 
% [ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(biasDiff(:,1));%phantom low vs phantom control low
% [ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(biasDiff(:,2));%phantom med vs phantom control med
% [ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(biasDiff(:,3));%phantom high vs phantom control high
