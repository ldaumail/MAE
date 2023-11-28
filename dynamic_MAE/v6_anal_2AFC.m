names = {'sub-01'}; %1237910
 
name = names{1};
version = 'v6';
%if contains(name, s_f_names) %check if subject had the slow condition first or the fast one
folderDir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/data',...
        sprintf('/%s/%s_%s_2AFC',version,name,version), '/');
    
fileNames = dir(folderDir);
fileNames = fileNames(contains({fileNames.name},'MAE'));
fileDir = fileNames.folder;
exptDat = load(strcat(fileNames(1).folder,'/',fileNames(1).name));
exptDat = exptDat.ex;
conditions = exptDat.conds;
condNums = exptDat.condShuffle;

 
teststarts = exptDat.flipTime(exptDat.blockLength*exptDat.flipsPerSec,1:end)-exptDat.startRun; %exptDat.blockLength*exptDat.flipsPerSec%trial test start time relative to experiment start time
testends = exptDat.flipTime(end,1:end)-exptDat.startRun+exptDat.ITI1;  %trial test end time relative to experiment start time
responseType = [1, 2];

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
respFreq = squeeze(sum(respDat,2)); %n = response number (only one response per row across all response types), resps(n) = response type,


%% Look at frequency of the responses for each condition

 yvar = 100*[respFreq(:,1)./sum(respFreq(:,1)), respFreq(:,2)./sum(respFreq(:,2)), respFreq(:,5)./sum(respFreq(:,5)), respFreq(:,6)./sum(respFreq(:,6)), respFreq(:,3)./sum(respFreq(:,3)), respFreq(:,4)./sum(respFreq(:,4))];
%yvar =[respFreq(:,1), respFreq(:,2), respFreq(:,5), respFreq(:,6), respFreq(:,3), respFreq(:,4)];
percepts =  {'Probe going up','Probe going down'};

avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Dim bars match bg','(Phantom 1) up'), sprintf('%s\\newline%s\\newline%s\n','Dim bars match bg','(Phantom 1) down'), ...
    sprintf('%s\\newline%s\\newline%s\n','Light bars match bg',' lum (Phantom 2) up'),sprintf('%s\\newline%s\\newline%s\n','Light bars match bg',' lum (Phantom 2) down'), ...
   sprintf('%s\\newline%s\\newline%s\n','Inducer bar lum does not',' match bg (Control) up'),sprintf('%s\\newline%s\\newline%s\n','Inducer bar lum does not',' match bg (Control) down'),... % sprintf('%s\\newline%s\\newline%s\n','Inducer bar lum does not',' match bg single (Control)'),
    'No inducer'}; 
ylab = {'Proportion of occurence of each percept (%)'};
ylims = [0 100];
barPlot(yvar,avgConditions, percepts, ylab, ylims)
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('proportion_percepts_across_conditions_%s_%s_sp6_test4_2AFC.png', name, version)));
