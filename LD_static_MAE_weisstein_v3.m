%MAE Test stimulus with constant contrast v11

names = {'sub-01', 'sub-02', 'sub-03', 'sub-04','sub-05','sub-06','sub-07'}; %
%'sub-05' and 'sub-06' did not experience any MAE effect

version = 'v3';
responseType = [1 2 3]; %1 = Rightward motion, 2= Leftward motion, 3= No motion
resptimes = nan(4, 4, length(names)); %nan(max(exptDat.repsPerRun), length(exptDat.repsPerRun), length(names)); %response time data better organized: 1: trial number, 2: condition number, 3: subjects number
respDir = nan(4, 4, length(names)); 

for i =1:length(names)
    name = names{i};
   
    %if contains(name, s_f_names) %check if subject had the slow condition first or the fast one
    folderDir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/data/static_MAE',...
        sprintf('/%s/%s_%s',version,name,version), '/');
    
    fileNames = dir(folderDir);
    fileNames = fileNames(contains({fileNames.name},'MAE'));
    fileDir = fileNames.folder;
    exptDat = load(strcat(fileNames(1).folder,'/',fileNames(1).name));
    exptDat = exptDat.ex;
    conditions = exptDat.conds;
    condNums = exptDat.condShuffle;
    
    teststarts = exptDat.flipTime(exptDat.blockLength*exptDat.flipsPerSec,1:end)-exptDat.startRun; %exptDat.blockLength*exptDat.flipsPerSec%trial test start time relative to experiment start time
    testends = exptDat.responseTimes;% - (exptDat.flipTime(1,1:39)-exptDat.startRun);
    
    
    cnt = zeros(length(unique(condNums)),1);
    
    for c =1:length(condNums)
        condNum = condNums(c);
        cnt(condNum) = cnt(condNum)+1; %counting rep number
        tstart = teststarts(c);
        % Initialize tend to NaN
        
        tend = NaN;
        resp = NaN;
        if c < length(condNums)
            resp = exptDat.resp((testends > tstart & testends < teststarts(c + 1)));
            resp = resp(1);
            for t = 1:length(testends)
                
                if (testends(t) > tstart && testends(t) < teststarts(c + 1))
                    tend = testends(t);
                    
                    break
                end
                
            end
        elseif c == length(condNums)
            resp = exptDat.resp(testends > tstart);
            resp = resp(1);
            for t = 1:length(testends)
                if testends(t) > tstart
                    tend = testends(t);
                    
                    break
                end
            end
        end
        % Assign tend to respDat, or set to NaN if not found
        if ~isnan(tend)
            resptimes(cnt(condNum), condNum,i) = tend-tstart;
            respDir(cnt(condNum), condNum,i) = resp;
        else
            resptimes(cnt(condNum), condNum,i) = 0;
            respDir(cnt(condNum), condNum,i) = 0;
        end
    end
    
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

%% Parse out same/different directions

same_diff = nan(4,length(conditions), length(names));
for i = 1:length(names)
    for c =1:length(conditions)
        for t = 1:size(respDir,1)
                %same direction
            if contains(conditions(c), 'Right') && respDir(t,c,i) == 1
                same_diff(t,c,i) = NaN;
            elseif contains(conditions(c), 'Left') && respDir(t,c,i) == 2
                same_diff(t,c,i) = NaN;
                %ambiguous direction/no MAE
            elseif respDir(t,c,i) == 3
                same_diff(t,c,i) = NaN;%0.5;
                %opposite direction
            elseif contains(conditions(c), 'Right') && respDir(t,c,i) == 2
                same_diff(t,c,i) = 1;
            elseif  contains(conditions(c), 'Left') && respDir(t,c,i) == 1
                same_diff(t,c,i) = 1;
            end
        end
    end
end



%% Look at "Phantom" vs "phantomcontrol" durations of MAE
resptimes = resptimes.*same_diff;
yvar = [[resptimes(:,1,:); resptimes(:,2,:)], [resptimes(:,3,:); resptimes(:,4,:)]]; 
% ydat = permute(ydat, [1,3,2]);
yvar = squeeze(nanmean(yvar,1))';
avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Phantom'),  ...
    sprintf('%s\\newline%s\\newline%s\n','Phantom','Control')...
   };
percepts = {''};

ylab = {'MAE duration (s)'};
ylims = [0 5];

groupsDotBarPlotSEM(yvar,avgConditions, percepts, ylab)
 %singleBarDotPlotSEM3(yvar,avgConditions, percepts, ylab, ylims)
%barPlotSD2(yvar,avgConditions, ylab, ylims)
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('MAE_duration_%s_weisstein.png',version)));

%% Percentage of trials in opposite direction vs not.

oppositeDir = nan(4,length(conditions), length(names));
sameDir = nan(4,length(conditions), length(names));
ambigDir = nan(4,length(conditions), length(names));
for i = 1:length(names)
    for c =1:length(conditions)
        for t = 1:size(respDir,1)
                %same direction
            if contains(conditions(c), 'Right') && respDir(t,c,i) == 1
                sameDir(t,c,i) = 1;
            elseif contains(conditions(c), 'Left') && respDir(t,c,i) == 2
                sameDir(t,c,i) = 1;
                %ambiguous direction/no MAE
            elseif respDir(t,c,i) == 3
                ambigDir(t,c,i) = 1;
                %opposite direction
            elseif contains(conditions(c), 'Right') && respDir(t,c,i) == 2
                oppositeDir(t,c,i) = 1;
            elseif  contains(conditions(c), 'Left') && respDir(t,c,i) == 1
                oppositeDir(t,c,i) = 1;

            end
        end
    end
end

combOppositeDir =[[oppositeDir(:,1,:);oppositeDir(:,2,:)],[oppositeDir(:,3,:);oppositeDir(:,4,:)]];
combAmbigDir =[[ambigDir(:,1,:);ambigDir(:,2,:)],[ambigDir(:,3,:);ambigDir(:,4,:)]];
combSameDir =[[sameDir(:,1,:);sameDir(:,2,:)],[sameDir(:,3,:);sameDir(:,4,:)]];

percOppDir = nan(length(names),2);
percSameDir = nan(length(names),2);
percAmbigDir = nan(length(names),2);


for i =1:length(names)
    for c =1:2
        percOppDir(i,c) = 100*sum(~isnan(combOppositeDir(:,c,i)))/length(combOppositeDir(:,c,i));
        percAmbigDir(i,c) = 100*sum(~isnan(combAmbigDir(:,c,i)))/length(combAmbigDir(:,c,i));
        percSameDir(i,c) = 100*sum(~isnan(combSameDir(:,c,i)))/length(combSameDir(:,c,i));
    end
end

conds = {'Phantom','Phantom control'};
respType = {'Opposite', 'Ambiguous', 'Same'};
yval = cat(3, percOppDir, percAmbigDir,percSameDir);
ylab = 'Percent of occurence (%)';
% groupsBarPlotSEM(permute(yval, [3 1 2]),conds, respType, ylab) 
singleBarLinePlotSEM3D(yval,conds,respType, ylab, [0 100])
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('MAE_percent_type_%s_weisstein_subjs.png',version)));

%% Plot histogram of Bias score

%bias
score=[0; 0.5; 1]; %score each percept type: 0 = same, 0.5 = ambiguous, 1 = opposite

bias = squeeze(sum(combSameDir*score(1),1,'omitnan') + sum(combAmbigDir*score(2),1,'omitnan') + sum(combOppositeDir*score(3),1,'omitnan'));

totalResps = squeeze(sum(combSameDir,1,'omitnan') + sum(combAmbigDir,1,'omitnan') + sum(combOppositeDir,1,'omitnan'));

percentBias = 100*bias./totalResps;

yvar = percentBias;
conds = {'Phantom','Phantom control'};
xlab = {'Percent bias (%)'};
ylab = {'Frequency'};
percentHistogram(yvar,conds, xlab, ylab,[0 3])
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('sMAE_percent_bias_hist_%s.png',version)));

% figure();
% subplot(2,1,1)
% histogram(percentBias(1,:),'BinWidth',5)
% subplot(2,1,2)
% histogram(percentBias(2,:),'BinWidth',5)



