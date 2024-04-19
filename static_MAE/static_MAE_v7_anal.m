names = {'sub-01','sub-02','sub-03', 'sub-04', 'sub-05', 'sub-06', 'sub-07', 'sub-08','sub-09', 'sub-10', 'sub-11', 'sub-12', 'sub-13', 'sub-14'}; %'sub-03'
version = 'v7';
responseType = [1, 2, 3];
respFreq = nan(length(responseType), 18, length(names));
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
        resps
        if ~isempty(resps)
            
            respDat(resps(end), cnt(condNum), condNum) = respDat(resps(end), cnt(condNum), condNum)+1; %resps(end) to take the last response in case the participant corrected a mistake
            
        end
    end
    respFreq(:,1:length(conditions),i) = squeeze(sum(respDat,2));
end

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
score=[0; 0.5;1]; %score each percept type: 0 = same, 0.5 = ambiguous, 1 = opposite
bias = same_diff.*score;
percentBias = nan(size(bias,2), length(names));
for i = 1:length(names)
    for c =1:size(bias,2)
        
        percentBias(c,i) = 100*sum(bias(:,c,i))/sum(same_diff(:,c,i));
        %     fprintf('Percent Bias %d\n', percentBias(c));
    end
end

contLevel = {'Low', 'Medium', 'High'};
conds = {'Full Grating','Phantom','Phantom Control'};

% condPercentBias = [percentBias(1,:); percentBias(2,:); percentBias(3,:); percentBias(4,:); percentBias(6,:); percentBias(8,:);  percentBias(5,:); percentBias(7,:); percentBias(9,:)];
condPercentBias = [percentBias(3,:); percentBias(6,:); percentBias(9,:); percentBias(1,:); percentBias(4,:); percentBias(7,:);  percentBias(2,:); percentBias(5,:); percentBias(8,:)];
orderedConds = {'Full Low', 'Full Med', 'Full High', 'Phantom Low', 'Phantom Med', 'Phantom High',  'Phantom Control Low', 'Phantom Control Med', 'Phantom Control High'};

yvar = reshape(condPercentBias, length(contLevel),length(conds),length(names));

yval = permute(yvar, [3 1 2]); %subj resp x contrast level x cond

ylab = {'Percent bias (%)'};
ylims = [40 110];

%  singleBarLinePlotSEM(condPercentBias,orderedConds, ylab, ylims)%  singleBarPlot(yvar(:,1)',avgConditions, {'bias'}, ylab, ylims);
MAEBarPlotSEM(yval,conds, ylab, ylims, contLevel)

% singleBarDotPlotSEM3(yvar,conds, {'bias'}, ylab, ylims);
 plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('staticMAE_percent_bias_%s.png', version)));

%%

