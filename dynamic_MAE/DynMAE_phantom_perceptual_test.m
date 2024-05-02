%% Analyze perceptual report
%% Loic Daumail - 03/21/2024
version = 'v13';
names = {'sub-01','sub-02','sub-03', 'sub-04', 'sub-05', 'sub-06', 'sub-07','sub-08','sub-09', 'sub-10', 'sub-11', 'sub-12', 'sub-13', 'sub-14', 'sub-15'};%

resps = nan(2,12,length(names));
for i =1:length(names)
    name = names{i};
    folderDir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/data',...
        sprintf('/dynamic_MAE/%s/%s_%s',version,name,version), '/');
    fileNames = dir(folderDir);
    fileNames = string({fileNames(contains({fileNames.name},'percept')).name});

    exptName = strcat(folderDir,fileNames(contains(fileNames,'percept')));
    exptDat = load(exptName);
    exptDat = exptDat.ex;
    conditions = exptDat.conds;
    condNums = exptDat.condShuffle;
    
    trConds = exptDat.condResp; %time at which a cycle of a condition starts. May be multiple cycles per trial for allowing time for participant to respond
    trChange = [1, diff(trConds)]; %Focus on cycles for which the condition changed = trial starts
    cycleStarts = exptDat.trialStartTime-exptDat.startRun; %Time that a cycle starts %exptDat.flipTime(exptDat.blockLength*exptDat.flipsPerSec,1:end)-exptDat.startRun; %exptDat.blockLength*exptDat.flipsPerSec%trial test start time relative to experiment start time
    trstarts = cycleStarts(cycleStarts & trChange); %Get trial start time
    trRespsT = exptDat.responseTimes; % Get response time
    trResps = exptDat.resp; %Get response type
    
    cnt = zeros(length(unique(condNums)),1);
    c = 0;
    for t =1:length(trstarts)
        condNum = condNums(t);
        cnt(condNum) = cnt(condNum)+1; %counting rep number
        if t < length(trstarts)
            if nnz((trRespsT > trstarts(t) & trRespsT < trstarts(t + 1)))
                resp = exptDat.resp((trRespsT > trstarts(t) & trRespsT < trstarts(t + 1)));
                resp = resp(end);
                resps(cnt(condNum),condNum,i) = resp;
            end
        elseif t == length(trstarts) & nnz((trRespsT > trstarts(t)))
                resp = exptDat.resp((trRespsT > trstarts(t)));
                resp = resp(end);
                resps(cnt(condNum),condNum,i) = resp;
        end
    end
    
end


yvar = squeeze(nanmean([[resps(:,1,:);resps(:,2,:)],[resps(:,5,:);resps(:,6,:)],[resps(:,9,:);resps(:,10,:)],[resps(:,3,:);resps(:,4,:)],[resps(:,7,:);resps(:,8,:)],[resps(:,11,:);resps(:,12,:)]],1))';

%% Plot
yval = reshape(yvar,size(yvar,1),3, 2); % %subj resp x contrast level x cond 

avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Low Contrast','Phantom'),  ...
    sprintf('%s\\newline%s\\newline%s\n','Medium Contrast','Phantom')...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast','Phantom')...
    sprintf('%s\\newline%s\\newline%s\n','Low Contrast', 'Phantom Control'),  ...
    sprintf('%s\\newline%s\\newline%s\n','Medium Contrast','Phantom Control')...
    sprintf('%s\\newline%s\\newline%s\n','High Contrast','Phantom Control')...
   };

percepts = {''};

ylab = {'Percept strength (avg score)'};
ylims = [0 5];
leg = {'Low', 'Medium', 'High'};
% singleBarLinePlotSEM(yvar',avgConditions, ylab, ylims)
MAEBarPlotSEM(yval,{'Phantom','Phantom Control'}, ylab, ylims, leg)
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('perceptual_report_dynamicMAE.svg')));


%% Stats

phcond = {'Phantom','Phantom','Phantom','PhantomControl','PhantomControl','PhantomControl'};
contLevels = {'Low';'Med';'High';'Low';'Med';'High'};

contrasts = [];
phantoms = [];
ph = [1 1 1 2 2 2]';
cont =[1 2 3 1 2 3]';

for i =1:length(contLevels)
    contrasts = [contrasts; repmat(cont(i),length(names),1)];
    phantoms = [phantoms; repmat(ph(i),length(names),1)];
end

subjectsIdx = repmat(names',length(contLevels),1);%repmat((1:length(names))',length(condNames),1);

data = reshape(yvar, [size(yvar,1)*size(yvar,2),1]);
tbl = table(subjectsIdx, data, contrasts, phantoms,'VariableNames',{'SubjectIndex','Response','Contrast','Phantom'});
lme = fitlme(tbl,'Response~Contrast*Phantom+(1|SubjectIndex)+(Contrast-1|SubjectIndex)+(Phantom-1|SubjectIndex)'); %

[pVal, F, R] = coefTest(lme);

%% ttests

[ttestMean(1), Pval(1),~,Stats(1).stats] = ttest(yvar(:,1),yvar(:,2));%phantom low vs phantom med
[ttestMean(2), Pval(2),~,Stats(2).stats] = ttest(yvar(:,1),yvar(:,3));%phantom low vs phantom high
[ttestMean(3), Pval(3),~,Stats(3).stats] = ttest(yvar(:,2),yvar(:,3));%phantom med vs phantom high

[ttestMean(4), Pval(4),~,Stats(4).stats] = ttest(yvar(:,4),yvar(:,5));%phantom control low vs phantom control med
[ttestMean(5), Pval(5),~,Stats(5).stats] = ttest(yvar(:,4),yvar(:,6));%phantom control low vs phantom control high
[ttestMean(6), Pval(6),~,Stats(6).stats] = ttest(yvar(:,5),yvar(:,6));%phantom control med vs phantom control high

