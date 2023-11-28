names = {'sub-06'}; %1237910
 
name = names{1};
version = 'v9';
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
responseType = [1, 2, 3];

respDat = zeros(length(responseType), max(exptDat.repsPerRun), exptDat.numConds); %response time data better organized
% respFreq = nan(length(responseType), max(exptDat.repsPerRun), exptDat.numConds); %response frequency data better organized

cnt = zeros(length(condNums),1);

for i =1:length(condNums)
    condNum = condNums(i);
    cnt(condNum) = cnt(condNum)+1; %counting rep number
    tstart = teststarts(i);
    tend = testends(i);
    respTimes = exptDat.responseTimes(exptDat.responseTimes > tstart &  exptDat.responseTimes < tend);
    resps = exptDat.correctResp(exptDat.responseTimes > tstart &  exptDat.responseTimes < tend);
    resps
    if ~isempty(resps)
        
        respDat(resps(end), cnt(condNum), condNum) = respDat(resps(end), cnt(condNum), condNum)+1; %resps(end) to take the last response in case the participant corrected a mistake
        
    end
end
respFreq = squeeze(sum(respDat,2)); 

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

same_diff = zeros(3,length(conditions)/2);
for c =1:length(conditions)
    for t = 1:size(respFreq,1)
        %same direction
        if contains(conditions(c), 'Up') && t == 1
            condcnt = ceil(c/2);
            same_diff(1,condcnt) = same_diff(1,condcnt)+ respFreq(t,c);
        elseif contains(conditions(c), 'Down') && t == 3
            condcnt = ceil(c/2);
            same_diff(1,condcnt) = same_diff(1,condcnt)+ respFreq(t,c);
        %ambiguous direction
        elseif t == 2 
            condcnt = ceil(c/2);
            same_diff(2,condcnt) = same_diff(2,condcnt)+ respFreq(t,c);
        %opposite direction
        elseif contains(conditions(c), 'Up') && t == 3
            condcnt = ceil(c/2);
            same_diff(3,condcnt) = same_diff(3,condcnt)+ respFreq(t,c);
        elseif  contains(conditions(c), 'Down') && t == 1
            condcnt = ceil(c/2);
            same_diff(3,condcnt) = same_diff(3,condcnt)+ respFreq(t,c);
        end   
    end
end


yvar = 100*[same_diff(:,1)./sum(same_diff(:,1)), same_diff(:,3)./sum(same_diff(:,3)), same_diff(:,5)./sum(same_diff(:,5)),...
     same_diff(:,2)./sum(same_diff(:,2)), same_diff(:,4)./sum(same_diff(:,4)), same_diff(:,6)./sum(same_diff(:,6)),...
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
barPlot2(yvar,avgConditions, percepts, ylab, ylims)
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_same_opposite_%s_%s_sp4_test2point67.png', name, version)));
