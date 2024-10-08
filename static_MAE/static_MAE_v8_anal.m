%Loic Daumail


names = {'sub-01','sub-02','sub-03', 'sub-04','sub-05'};%
version = 'v8';
responseType = [1, 2, 3];
respFreq = nan(length(responseType), 6, length(names));
respDur = nan(length(responseType),3, 6, length(names));
% gapSize = [];
for i =1:length(names)
    name = names{i};
   
    folderDir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/data',...
        sprintf('/static_MAE/%s/%s_%s',version,name,version), '/');
    
    fileNames = dir(folderDir);
    fileNames = fileNames(contains({fileNames.name},sprintf('MAE_%s',version)));
    fileDir = fileNames.folder;
    exptDat = load(strcat(fileNames(1).folder,'/',fileNames(1).name));
    exptDat = exptDat.ex;
    conditions = exptDat.conds;
    condNums = exptDat.condShuffle;
    
    teststarts = exptDat.flipTime(exptDat.blockLength*exptDat.flipsPerSec,1:end)-exptDat.startRun; %exptDat.blockLength*exptDat.flipsPerSec%trial test start time relative to experiment start time
    
    releaseTimes = exptDat.releaseTimes(exptDat.releaseTimes>0); 
    tends = [];
    for c =1:length(condNums)
        if c == length(condNums) && teststarts(c)< releaseTimes(end)
            clear tend
            tend = releaseTimes(end);
            tends = [tends, tend];
        elseif c <= length(condNums)
            clear tend
            tend = releaseTimes(teststarts(c)< releaseTimes & releaseTimes<teststarts(c+1));
            if ~isempty(tend)
                tends = [tends, tend];
            elseif isempty(tend)
                tends = [tends, NaN];
            end
        end
        
    end

    respDat = zeros(length(responseType), max(exptDat.repsPerRun), exptDat.numConds); %response time data better organized
    
    cnt = zeros(exptDat.numConds,1);
    
    for c =1:length(condNums)
        condNum = condNums(c);
        cnt(condNum) = cnt(condNum)+1; %counting rep number
        tstart = teststarts(c);
        tend = tends(c);

        respTimes = exptDat.pressTimes(exptDat.pressTimes > tstart &  exptDat.pressTimes < tend);
        resps = exptDat.resp(exptDat.pressTimes > tstart &  exptDat.pressTimes < tend);
%         resps
        if ~isempty(resps) 
            %Count response types
            respDat(resps(end), cnt(condNum), condNum) = respDat(resps(end), cnt(condNum), condNum)+1; %resps(end) to take the last response in case the participant corrected a mistake
            %Record Response duration from test onset time
            respDur(resps(end), cnt(condNum), condNum, i) = tend - tstart;    
        end
    end
    respFreq(:,1:length(conditions),i) = squeeze(sum(respDat,2));
end


%% Look at "same" vs "ambiguous" vs "different" direction percepts
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

%% plot : subjects scatter plots + bar dots grouped by percept
percepts =  flip({'Same','None','Opposite'});

ylab = {'Proportion of occurence of each percept (%)'};

contLevel = {'Low', 'Medium', 'High'};
conds = {'Full Grating','Phantom','Phantom Control'};
% condMAEDir is subjs * all conds * respType here
condMAEDir = [yvar(:,3,:), yvar(:,1,:), yvar(:,2,:)];
% orderedConds = {'Full Low', 'Phantom Low',  'Phantom Control Low'};

%Exclude subjects with less than 60 percent MAE
%avgMAE = nanmean(condMAEDir(:,1:3,3),2) < 60;


yvar2 = flip(condMAEDir,3);
ylims = flip([0 105; 0 105; 0 105]);
lowCtMAEBarPlotSEM2(yvar2,conds, percepts, ylab,ylims)
 plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
 saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_same_opposite_staticMAE_%s_group_subplot.png', version)));

 
 %%
 saveDir = "/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/curated_data";
mkdir(saveDir);
save(strcat(saveDir,"/static_MAE_dir.mat"), "condMAEDir")
 
%%
%Line plots

percepts =  flip({'Same','None','Opposite'});
ylab = {'Proportion of occurence of each percept (%)'};
% contLevel = {'Low', 'Medium', 'High'};
% conds = {'Full Grating','Phantom','Phantom Control'};
% condMAEDir is subjs * all conds * respType here
condMAEDir = [yvar(:,3,:), yvar(:,1,:), yvar(:,2,:)];
orderedConds = {'Full Grating','Phantom','Phantom Control'};


%yvar2 = flip(reshape(condMAEDir, length(names), length(contLevel), length(conds), length(responseType)),4);
ylims = flip([0 105; 0 105; 0 105]);

%same
singleBarLinePlotSEM(condMAEDir(:,:,1)',orderedConds, ylab, ylims(1,:))
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');

 saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_same_staticMAE_%s_lineplot.png', version)));

 %opposite
 singleBarLinePlotSEM(condMAEDir(:,:,3)',orderedConds, ylab, ylims(1,:))
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');

 saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_opposite_staticMAE_%s_lineplot.png', version)));

 %% stats

yvar = condMAEDir;

% Same 

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1,1),yvar(:,2,1));%Full Low vs Full medium
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,1,1),yvar(:,3,1));%Full Low vs Full High
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,2,1),yvar(:,3,1));%Full Medium vs Full high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4,1),yvar(:,5,1));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,4,1),yvar(:,6,1));%phantom low vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,5,1),yvar(:,6,1));%phantom med vs phantom high

[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,7,1),yvar(:,8,1));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,7,1),yvar(:,9,1));%phantom control low vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,8,1),yvar(:,9,1));%phantom control med vs phantom control high


% None
[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1,2),yvar(:,2,2));%Full Low vs Full medium
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,1,2),yvar(:,3,2));%Full Low vs Full High
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,2,2),yvar(:,3,2));%Full Medium vs Full high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4,2),yvar(:,5,2));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,4,2),yvar(:,6,2));%phantom low vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,5,2),yvar(:,6,2));%phantom med vs phantom high

[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,7,2),yvar(:,8,2));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,7,2),yvar(:,9,2));%phantom control low vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,8,2),yvar(:,9,2));%phantom control med vs phantom control high

% Opposite
[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1,3),yvar(:,2,3));%Full Low vs Full medium
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,1,3),yvar(:,3,3));%Full Low vs Full High
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,2,3),yvar(:,3,3));%Full Medium vs Full high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4,3),yvar(:,5,3));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,4,3),yvar(:,6,3));%phantom low vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,5,3),yvar(:,6,3));%phantom med vs phantom high

[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,7,3),yvar(:,8,3));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,7,3),yvar(:,9,3));%phantom control low vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,8,3),yvar(:,9,3));%phantom control med vs phantom control high

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

[pVal, F, DF1, DF2] = coefTest(lme);

% 
%% No MAE
data = reshape(condMAEDir(:,:,2), [size(condMAEDir(:,:,2),1)*size(condMAEDir(:,:,2),2),1]);
tbl = table(subjectsIdx', data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, DF1, DF2] = coefTest(lme);
% 
% %% opposite
data = reshape(condMAEDir(:,:,3), [size(condMAEDir(:,:,3),1)*size(condMAEDir(:,:,3),2),1]);
tbl = table(subjectsIdx', data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, DF1, DF2] = coefTest(lme);
%% MAE Bias

% respFreq = squeeze(sum(respDat,2));

same_diff = zeros(length(responseType),length(conditions)/2, length(names));
for i = 1:length(names)
    for c =1:length(conditions)
        for t = 1:size(respFreq,1) %response type
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
score=[-1; 0;1]; %score each percept type: 0 = same, 0.5 = ambiguous, 1 = opposite
bias = same_diff.*score;
percentBias = nan(length(names),size(bias,2));
for i = 1:length(names)
    for c =1:size(bias,2)
        
        percentBias(i,c) = 100*sum(bias(:,c,i))/sum(same_diff(:,c,i));
        %     fprintf('Percent Bias %d\n', percentBias(c));
    end
end

contLevel = {'Low'};
conds = {'Full Grating','Phantom','Phantom Control'};

condPercentBias = [percentBias(:,3), percentBias(:,1), percentBias(:,2)];

yvar = condPercentBias;


ylab = {'Percent bias (%)'};
ylims = [0 100];

%  singleBarLinePlotSEM(condPercentBias',orderedConds, ylab, [-50 100])%  singleBarPlot(yvar(:,1)',avgConditions, {'bias'}, ylab, ylims);

MAEctBarPlotSEM(yvar,conds, ylab, ylims,contLevel)

% singleBarDotPlotSEM3(yvar,conds, {'bias'}, ylab, ylims);
 plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('staticMAE_percent_bias_subtract_%s.png', version)));

%% Stats
saveDir = "/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/curated_data";
mkdir(saveDir);
save(strcat(saveDir,"/static_MAE_percent_bias.mat"), "condPercentBias")


%% Inducer type vs contrast interaction

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

[pVal, F, DF1, DF2] = coefTest(lme);


%% ttests
yvar = condPercentBias;

% [ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,4),yvar(:,7));%phantom low vs phantom control low
% [ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,5),yvar(:,8));%phantom med vs phantom control med
% [ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,6),yvar(:,9));%phantom high vs phantom control high
% [ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4),yvar(:,5));%phantom low vs phantom med
% [ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,5),yvar(:,6));%phantom med vs phantom high
% [ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,4),yvar(:,6));%phantom low vs phantom high
% [ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,7),yvar(:,8));%phantom control low vs phantom control med
% [ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,8),yvar(:,9));%phantom control med vs phantom control high
% [ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,7),yvar(:,9));%phantom control low vs phantom control high

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1),yvar(:,2));%Full Low vs Full medium
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,1),yvar(:,3));%Full Low vs Full High
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,2),yvar(:,3));%Full Medium vs Full high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4),yvar(:,5));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,4),yvar(:,6));%phantom low vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,5),yvar(:,6));%phantom med vs phantom high

[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,7),yvar(:,8));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,7),yvar(:,9));%phantom control low vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,8),yvar(:,9));%phantom control med vs phantom control high

% 
%% Response duration

respCondDur(:,:,1,:) = [respDur(:,:,1,:), respDur(:,:,2,:)];
respCondDur(:,:,2,:) = [respDur(:,:,3,:), respDur(:,:,4,:)];
respCondDur(:,:,3,:) = [respDur(:,:,5,:), respDur(:,:,6,:)];
respCondDur(:,:,4,:) = [respDur(:,:,7,:), respDur(:,:,8,:)];
respCondDur(:,:,5,:) = [respDur(:,:,9,:), respDur(:,:,10,:)];
respCondDur(:,:,6,:) = [respDur(:,:,11,:), respDur(:,:,12,:)];
respCondDur(:,:,7,:) = [respDur(:,:,13,:), respDur(:,:,14,:)];
respCondDur(:,:,8,:) = [respDur(:,:,15,:), respDur(:,:,16,:)];
respCondDur(:,:,9,:) = [respDur(:,:,17,:), respDur(:,:,18,:)];


same_diff = zeros(length(responseType),max(exptDat.repsPerRun)*2, length(conditions)/2, length(names));
for i = 1:length(names)
    for c =1:length(conditions)/2
        for r = 1:max(exptDat.repsPerRun)*2
            for t = 1:length(responseType) %response type % Right arrow was 1, Left arrow was 2 and Down arrow was 3
                %same direction
                if contains(conditions(c), 'Right') && t ==1 && ~isnan(respCondDur(t,r,c,i))
                    same_diff(t,r,c,i) = 0;
                elseif contains(conditions(c), 'Left') && t == 2 && ~isnan(respCondDur(t,r,c,i))

                    same_diff(t,r,c,i) = 0;
                    %No motion aftereffect
                elseif t == 3 && ~isempty(respCondDur(t,r,c,i))
                    same_diff(t,r,c,i) = 0;
                    %opposite direction
                elseif contains(conditions(c), 'Right') && t == 2 && ~isnan(respCondDur(t,r,c,i))
                    same_diff(t,r,c,i) = 1;
                elseif  contains(conditions(c), 'Left') && t == 1 && ~isnan(respCondDur(t,r,c,i))
                    same_diff(t,r,c,i) = 1;
                    
                end
                
            end
        end
    end
end


selectRespDurs = respCondDur.*same_diff;

lineRespDur = squeeze(nanmean(selectRespDurs,1)); %Since there should be one value per column (one response per trial), is allows to reduce the dimensionality of the matrix from all response types


linearAvgRespDurs = squeeze(nanmean(lineRespDur,1))';


contLevel = {'Low', 'Medium', 'High'};
conds = {'Full Grating','Phantom','Phantom Control'};

responseDurations = [linearAvgRespDurs(:,3), linearAvgRespDurs(:,6), linearAvgRespDurs(:,9), linearAvgRespDurs(:,1), linearAvgRespDurs(:,4), linearAvgRespDurs(:,7),  linearAvgRespDurs(:,2), linearAvgRespDurs(:,5), linearAvgRespDurs(:,8)];
orderedConds = {'Full Low', 'Full Med', 'Full High', 'Phantom Low', 'Phantom Med', 'Phantom High',  'Phantom Control Low', 'Phantom Control Med', 'Phantom Control High'};

yvar = reshape(responseDurations, length(names), length(contLevel),length(conds));


ylab = {'MAE duration (s)'};
ylims = [0 2.5];

%  singleBarLinePlotSEM(condPercentBias,orderedConds, ylab, ylims)%  singleBarPlot(yvar(:,1)',avgConditions, {'bias'}, ylab, ylims);
MAEBarPlotSEM(yvar,conds, ylab, ylims, contLevel)

% singleBarDotPlotSEM3(yvar,conds, {'bias'}, ylab, ylims);
 plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('staticMAE_MAE_duration_zeroaverage_%s.svg', version)));   


%% Stats
saveDir = "/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/curated_data";
mkdir(saveDir);
save(strcat(saveDir,"/static_MAE_duration.mat"), "responseDurations")
%% Inducer type vs contrast interaction

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


data = reshape(responseDurations, [size(responseDurations,1)*size(responseDurations,2),1]);
tbl = table(subjectsIdx, data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);
%% ttests
yvar = responseDurations;

% [ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,4),yvar(:,7));%phantom low vs phantom control low
% [ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,5),yvar(:,8));%phantom med vs phantom control med
% [ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,6),yvar(:,9));%phantom high vs phantom control high
% [ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4),yvar(:,5));%phantom low vs phantom med
% [ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,5),yvar(:,6));%phantom med vs phantom high
% [ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,4),yvar(:,6));%phantom low vs phantom high
% [ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,7),yvar(:,8));%phantom control low vs phantom control med
% [ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,8),yvar(:,9));%phantom control med vs phantom control high
% [ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,7),yvar(:,9));%phantom control low vs phantom control high
[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1),yvar(:,2));%Full Low vs Full medium
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,1),yvar(:,3));%Full Low vs Full High
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,2),yvar(:,3));%Full Medium vs Full high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4),yvar(:,5));%phantom low vs phantom med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,4),yvar(:,6));%phantom low vs phantom high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,5),yvar(:,6));%phantom med vs phantom high

[ttestMean(7), Pval(7),~,Stats(7).stats] = ttest(yvar(:,7),yvar(:,8));%phantom control low vs phantom control med
[ttestMean(8), Pval(8),~,Stats(8).stats] = ttest(yvar(:,7),yvar(:,9));%phantom control low vs phantom control high
[ttestMean(9), Pval(9),~,Stats(9).stats] = ttest(yvar(:,8),yvar(:,9));%phantom control med vs phantom control high

