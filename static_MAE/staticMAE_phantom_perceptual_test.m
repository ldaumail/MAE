version = 'v7';
names = {'sub-01','sub-02','sub-03', 'sub-04', 'sub-05', 'sub-06', 'sub-07', 'sub-08','sub-09'};%

resps = nan(4,3,length(names));
for i =1:length(names)
    name = names{i};
    folderDir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/data',...
        sprintf('/static_MAE/%s/%s_%s',version,name,version), '/');
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


avgResps = nanmean(resps,1);

%%

yvar = squeeze(avgResps); 
yval = [yvar(1,:),yvar(2,:); yvar(5,:),yvar(6,:); yvar(9,:),yvar(10,:); yvar(3,:),yvar(4,:); yvar(7,:),yvar(8,:); yvar(11,:),yvar(12,:)];

yval = reshape(yval,3, 2,size(yval,2)); % contrast level x cond x subj resps   
yval = permute(yval,[3 1 2]);%subj resp x contrast level x cond

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
% singleBarLinePlotSEM(yval,avgConditions, ylab, ylims)
MAEBarPlotSEM(yval,{'Phantom','Phantom Control'}, ylab, ylims, leg)
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('perceptual_report_staticMAE.png')));
