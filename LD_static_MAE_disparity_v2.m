names = {'sub-01', 'sub-02', 'sub-03', 'sub-04'}; %'sub-03'
version = 'v3';
responseType = [1 2 3]; %1 = Rightward motion, 2= Leftward motion, 3= No motion
resptimes = nan(4, 4, length(names)); %nan(max(exptDat.repsPerRun), length(exptDat.repsPerRun), length(names)); %response time data better organized: 1: trial number, 2: condition number, 3: subjects number
respDir = nan(4, 4, length(names)); 

for i =1:length(names)
    name = names{i};
   
    %if contains(name, s_f_names) %check if subject had the slow condition first or the fast one
    folderDir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/data/static_MAE',...
        sprintf('/stereo_%s/%s_%s',version,name,version), '/');
    
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
            for t = 1:length(testends)
                
                if (testends(t) > tstart && testends(t) < teststarts(c + 1))
                    tend = testends(t);
                    
                    break
                end
                
            end
        elseif c == length(condNums)
            resp = exptDat.resp(testends > tstart);
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
            respDir(cnt(condNum), condNum,i) = resp(1);
        else
            resptimes(cnt(condNum), condNum,i) = 0;
            respDir(cnt(condNum), condNum,i) = 0;
        end
    end
    
end

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

resptimes = resptimes.*same_diff;
yvar = [[resptimes(:,1,:); resptimes(:,2,:)], [resptimes(:,3,:); resptimes(:,4,:)]]; 
% ydat = permute(ydat, [1,3,2]);
% yvar = squeeze(nanmean(yvar,1))';
avgConditions = {sprintf('%s\\newline%s\\newline%s\n','Phantom', 'Front'),  ...
    sprintf('%s\\newline%s\\newline%s\n','Phantom','Back')...
   };
percepts = {''};

ylab = {'MAE duration (s)'};
ylims = [0 5];

%groupsDotBarPlotSEM(yvar,avgConditions, percepts, ylab)
% singleBarDotPlotSEM3(yvar,avgConditions, percepts, ylab, ylims)
barPlotSD2(yvar,avgConditions, ylab, ylims)
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('staticMAE_duration_%s_disparity_%s.png', name,version)));

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
respType = {'Opposite', 'Ambiguous'};%, 'Same'
yval = cat(3, percOppDir, percAmbigDir); %percSameDir
ylab = 'Percent of occurence (%)';
groupsBarPlotSEM(permute(yval, [1 3 2]),conds, respType, ylab)

plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/motion_after_effect/anal_plots/');
mkdir(plotdir);
saveas(gcf,strcat(plotdir, sprintf('MAE_percent_type_%s_disparity.png',version)));

