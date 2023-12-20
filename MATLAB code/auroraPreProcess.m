%{ 

auroraPreProcess.m

1. Initialization
* Dataset name, expected number of sensors, sensor measurement frequency
* Filter options: lowhz: butterworth filter low pass Hz

2. Clean raw data into usable format
* Import raw aurora output, fixes default timetable format, sets channel
names
* Save clean data

3. Interpolate and Filter
* Import clean data, interpolate missing data, filter with butterworth low
pass filter (plot channels for all three datasets)
* Save interpolated and filtered data

4. Inspection: Choose metric to calculate midpoint (inflection point)
* Compare the findchangepts results for midpoint (between A and B sets of
motion) using each metric (linear, std, mean, rms), for both interpolated +
filtered signal (using channel s1Tz)

5. Calculate and plot averaged midpoint
* Using the result from 4, graph the calculated midpoint on the filtered
channel data; Calculate the average midpoint using the channels s1Tx, s1Tz,
s1Ry, s1Rz.
* Review the averaged midpoint on the s1Tz graph

6. Separate front and back segments
* Clip the raw and filtered data at the midpoint (the data at the midpoint
being added to the front/'A' set)

7. Review motion sets for Error Percent
* Plots the raw and interpolated channels

8. Locate beginning of Dataset 'A'
* Plot local extrema for each channel in the dataset. Use the trough number
for the beginning of the motion in channels s1Tx and s1Tz.

9. Isolate Dataset 'A'
* Using the trough numbers from 8, average the location and split both the
raw and filtered 'front' data, into 'init' and 'A'. The location itself to
be contained in Dataset 'A'.

10. Review error in Dataset 'A'
__

Data Saved:
* CLEAN: Raw data in a usable format:
    (FEB 13 CLEAN DATA/aurora-clean-[experiment].csv)
* CLEAN [ERROR]: Raw data with a channel missing:
    (FEB 13 CLEAN DATA/aurora-clean-ERROR-[experiment].csv)
* INTERP: Interpolated data:
    (FEB 13 INTERP DATA/aurora-interp-feb-13-[experiment].csv)
* FILTERED: Filtered (butterworth LP) data:
    (FEB 13 FILTERED DATA/aurora-butter-feb-13-[experiment].csv)
* CLIPPED: Front/Back half of the raw & filtered data
    (FEB 13 CLIPPED DATA/aurora-[front/back]-[experiment].csv)
    (FEB 13 CLIPPED DATA/aurora-prefiltered-[front/back]-[experiment].csv
* Datasets: Init and dataset 'A', raw & filtered
    (Datasets/Raw: aurora-[experiment]-raw-[a/init].csv
    (Datasets: aurora-[experiment]-[a/init].csv

Nicole Devos for the WearME Lab, Western University
Version 2.3
July 20, 2023

%}

%% Initialization

close all
%clear all
%clc

% Data set name
% Available options for names:
% 
% test motions:
% "motor-trigger-test"; "cage-max-angle-test"; "side-side-motor-test";
% "side-side-test"; "side-side-trigger-test"; "zero-degree-motor-test";
% "zero-degree-test";
% 
% constant values:
% "base-x-z"; "base-y"; "cage-landmark"; "jt2-x"; "skull-landmark";
% "wheel-landmark";
% 
% Data sets missing sensor 2:
% "exp2-trigger-3"; "exp3-1"; "exp3-2"; "exp3-3"; "exp4-trigger-1";
% "exp4-trigger-2"; "exp4-trigger-3";

 name = 'exp1-1';
% name = 'exp1-2';
% name = 'exp1-3';
% name = 'exp2-trigger-1'; % Note: this one has a 43.1% error in the A set
% name = 'exp2-trigger-2';
% name = 'exp5-1';
% name = 'exp5-2';
% name = 'exp5-3';
% name = 'exp6-trigger-1';
% name = 'exp6-trigger-2';
% name = 'exp6-trigger-3';

numSensor = 2; % expected number of sensors for the data set
aurFreq = 40; % Aurora frequency 40 Hz

% Filter options
lowhz = 0.5; % Butterworth, 0.5 Hz

%% Import Raw, remove extra columns, save as "CLEAN" version
% Doesn't work in 2022a:
%{
aurVariableNames = ["s1Rz","s1Ry","s1Rx","s1Tx","s1Ty","s1Tz","s2Rz","s2Ry","s2Rx","s2Tx","s2Ty","s2Tz"];
altAurVariableNames = ["s1Rz","s1Ry","s1Rx","s1Tx","s1Ty","s1Tz"];

TTraw = [readtimetable(strcat('FEB 13 DATA/aurora-feb-13-',name,'.csv'),"SampleRate",aurFreq,"Range",'F:K') readtimetable(strcat('FEB 13 DATA/aurora-feb-13-',name,'.csv'),"SampleRate",aurFreq,"Range",'CT:CY')];
    TTraw.Time = duration(TTraw.Time,'Format','hh:mm:ss.SSSS');
    % Original data lists number of sensors
    numChannels = readvars(strcat('FEB 13 DATA/aurora-feb-13-',name,'.csv'),'Range','A1:A2');
    if numChannels == 1 % If one sensor
        TTraw.Properties.VariableNames = altAurVariableNames;
    else % If two sensors
        TTraw.Properties.VariableNames = aurVariableNames;
    end
    if numChannels ~= numSensor % If the dataset doesn't have the expected number of sensors
        disp(strcat('error in ',name,' dataset. Missing sensor data.'))
        %writetimetable(TTraw,strcat('FEB 13 CLEAN DATA/aurora-clean-feb-13-ERROR-',name,'.csv'));
    else
        %writetimetable(TTraw,strcat('FEB 13 CLEAN DATA/aurora-clean-feb-13-',name,'.csv'));
    end
%}

%% Interpolate, Filter, and plot - saved as "INTERP" and "FILTERED" data

% Raw data
TTclean = readtimetable(strcat('FEB 13 CLEAN DATA/aurora-clean-feb-13-',name,'.csv'));
plotAurora(TTclean,strcat("Individual measurements in the ",name," data"));

% Interpolated data
TTinterp = fixerrors(TTclean);
%writetimetable(TTinterp,strcat('FEB 13 INTERP DATA/aurora-interp-feb-13-',name,'.csv'));
plotAurora(TTinterp,strcat("Interpolated - Individual measurements in the ",name," data"));

% Filtered - butterworth
TTbutter = filterTT(TTinterp,aurFreq,lowhz);
%writetimetable(TTbutter,strcat('FEB 13 FILTERED DATA/aurora-butter-feb-13-',name,'.csv'));
plotAurora(TTbutter,strcat("Butterworth Filter @ 0.5 Hz - ",name));

%% Plotting findchangepts with FILTERED DATA
% Compare the raw, filtered, & raw + filtered changepoint using RMS, Mean,
% Linear, and Standard Deviation as the metrics.

% Testing std
pts.std = findchangepts(TTclean.s1Tz,'Statistic','std');
pts.stdF = findchangepts(TTbutter.s1Tz,'Statistic','std');

figure('Name','Raw vs Filtered Change points - std');
plot(TTclean.Time,TTclean.s1Tz)
title(strcat(name, ": Raw vs Filtered Change points - std"))
hold on
plot(TTbutter.Time,TTbutter.s1Tz,'r')
xline(TTclean.Time(pts.std),'-b',datestr(TTclean.Time(pts.std),'SS.FFF'));
xline(TTbutter.Time(pts.stdF),'--r',datestr(TTbutter.Time(pts.stdF),'SS.FFF'));

figure('Name','Raw Change Points - std');
findchangepts(TTclean.s1Tz,'Statistic','std');
figure('Name','Filtered Change Points - std');
findchangepts(TTbutter.s1Tz,'Statistic','std');

% Testing mean
pts.mean = findchangepts(TTclean.s1Tz,'Statistic','mean');
pts.meanF = findchangepts(TTbutter.s1Tz,'Statistic','mean');

figure('Name','Raw vs Filtered Change points - mean');
plot(TTclean.Time,TTclean.s1Tz)
title(strcat(name, ": Raw vs Filtered Change points - mean"))
hold on
plot(TTbutter.Time,TTbutter.s1Tz,'r')
xline(TTclean.Time(pts.mean),'-b',datestr(TTclean.Time(pts.mean),'SS.FFF'));
xline(TTbutter.Time(pts.meanF),'--r',datestr(TTbutter.Time(pts.meanF),'SS.FFF'));

figure('Name','Raw Change Points - mean');
findchangepts(TTclean.s1Tz,'Statistic','mean');
figure('Name','Filtered Change Points - mean');
findchangepts(TTbutter.s1Tz,'Statistic','mean');

% Testing rms
pts.rms = findchangepts(TTclean.s1Tz,'Statistic','rms');
pts.rmsF = findchangepts(TTbutter.s1Tz,'Statistic','rms');

figure('Name','Raw vs Filtered Change points - rms');
plot(TTclean.Time,TTclean.s1Tz)
title(strcat(name, ": Raw vs Filtered Change points - rms"))
hold on
plot(TTclean.Time,TTbutter.s1Tz,'r')
xline(TTclean.Time(pts.rms),'-b',datestr(TTclean.Time(pts.rms),'SS.FFF'));
xline(TTbutter.Time(pts.rmsF),'--r',datestr(TTbutter.Time(pts.rmsF),'SS.FFF'));

figure('Name','Raw Change Points - rms');
findchangepts(TTclean.s1Tz,'Statistic','rms');
figure('Name','Filtered Change Points - rms');
findchangepts(TTbutter.s1Tz,'Statistic','rms');

% Testing linear
pts.linear = findchangepts(TTclean.s1Tz,'Statistic','linear');
pts.linearF = findchangepts(TTbutter.s1Tz,'Statistic','linear');

figure('Name','Raw vs Filtered Change points - linear');
plot(TTclean.Time,TTclean.s1Tz)
title(strcat(name, ": Raw vs Filtered Change points - linear"))
hold on
plot(TTbutter.Time,TTbutter.s1Tz,'r')
xline(TTclean.Time(pts.linear),'-b',datestr(TTclean.Time(pts.linear),'SS.FFF'));
xline(TTbutter.Time(pts.linearF),'--r',datestr(TTbutter.Time(pts.linearF),'SS.FFF'));

figure('Name','Raw Change Points - linear');
findchangepts(TTclean.s1Tz,'Statistic','linear');
figure('Name','Filtered Change Points - linear');
findchangepts(TTbutter.s1Tz,'Statistic','linear');

%% Using the best findchangepoints metric; check the 12 channels
% Choose which is best in the previous section, then run it through the 12
% channels. Looking at s1Tx, s1Tz, s1Ry, s1Rz specifically.

% exp1-1: mean^           exp1-2: mean           exp1-3: mean
% exp2-trigger-1: mean*   exp2-trigger-2: mean 
% exp5-1: mean            exp5-2: mean            exp5-3: mean 
% exp6-trigger-1: mean    exp6-trigger-2: mean   exp6-trigger-3: mean
% * rms gets mess up when averaging
% ^ did  mean cos everything else was

 stat = 'mean';
% stat = 'std';
% stat = 'linear';
% stat = 'rms';

channel = {'s1Tx','s2Tx','s1Ty','s2Ty','s1Tz','s2Tz','s1Rx','s2Rx','s1Ry','s2Ry','s1Rz','s2Rz'};

figure('Name',strcat("Change points for ", name, " using '", stat, "'"));
sgtitle(strcat("Change points in dataset ", name, " using '", stat, "'"))

for i = 1:length(channel)
    pts.(channel{i})= findchangepts(TTbutter.(channel{i}),'Statistic',stat);
    subplot(6,2,i);
    plot(TTbutter.Time,TTbutter.(channel{i}))
    title(channel{i})
    hold on
    xline(TTbutter.Time(pts.(channel{i})));
end

% Check average of the above midpoints

% Calculate avg with the channels we are focusing one - using scientific
% rounding 
pts.mid = round(mean([pts.s1Tx,pts.s1Tz,pts.s1Ry,pts.s1Rz]),0,TieBreaker="even");

figure('Name',strcat("Averaged Midpoint for ", name, " using '", stat,"'"));
plot(TTclean.Time,TTclean.s1Tz)
title(strcat("Averaged midpoint for ",name," using '",stat,"'"))
hold on
plot(TTbutter.Time,TTbutter.s1Tz,'r')
xline(TTclean.Time(pts.mid),'-k',datestr(TTclean.Time(pts.mid),'SS.FFF'));

%% Clip and save at mid point; save as CLIPPED DATA
% Clippin includes the midpoint as within the first half of the motion
% (motion set 'A').

% Raw data
TTfront = TTclean(1:pts.mid,:);
%writetimetable(TTfront,strcat('FEB 13 CLIPPED DATA/aurora-front-feb-13-',name,'.csv'));
figure('Name',strcat("First half of motion for ", name, " (raw)"));
plot(TTfront.Time,TTfront.s1Tz)
title(strcat("First half of motion for ", name, " (raw)"))

TTback = TTclean((pts.mid+1):end,:);
%writetimetable(TTback,strcat('FEB 13 CLIPPED DATA/aurora-back-feb-13-',name,'.csv'));
figure('Name',strcat("Second half of motion for ", name, " (raw)"));
plot(TTback.Time,TTback.s1Tz)
title(strcat("Second half of motion for ", name, " (raw)"))

% Filtered data
TTbuttercF = TTbutter(1:pts.mid,:);
%writetimetable(TTbuttercF,strcat('FEB 13 CLIPPED DATA/aurora-prefiltered-front-feb-13-',name,'.csv'));
figure('Name',strcat("First half of motion for ", name, " (filtered)"));
plot(TTbuttercF.Time,TTbuttercF.s1Tz)
title(strcat("First half of motion for ", name, " (filtered)"))

TTbuttercB = TTbutter((pts.mid+1):end,:);
%writetimetable(TTbuttercB,strcat('FEB 13 CLIPPED DATA/aurora-prefiltered-back-feb-13-',name,'.csv'));
figure('Name',strcat("Second half of motion for ", name, " (filtered)"));
plot(TTbuttercB.Time,TTbuttercB.s1Tz)
title(strcat("Second half of motion for ", name, " (filtered)"))

% plot overlay
figure('Name',strcat("Split motion for ", name));
subplot(2,1,1)
plot(TTfront.Time,TTfront.s1Tz)
title(strcat("First half of motion for ", name))
% 875bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbgggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
hold on;
plot(TTbuttercF.Time,TTbuttercF.s1Tz)
legend('Raw','Filtered','Location','Southeast')

subplot(2,1,2)
plot(TTback.Time,TTback.s1Tz)
title(strcat("Second half of motion for ", name))
hold on;
plot(TTbuttercB.Time,TTbuttercB.s1Tz)
legend('Raw','Filtered','Location','Southeast')

%% Plot Aurora Channels + review percent error

% Plot raw channels (not interpolated)
plotAurora(TTfront,strcat("Front half of motion for ", name));
plotAurora(TTback,strcat("Back half of motion for ", name));

% Percent error
frontErr = checkError(TTfront);
disp(strcat("frontErr.s1Rx = ",num2str(frontErr.s1Rx*100),"%"));
if width(TTfront) == 12
    disp(strcat("frontErr.s2Rx = ",num2str(frontErr.s2Rx*100),"%"));
end
backErr = checkError(TTback);
disp(strcat("backErr.s1Rx = ",num2str(backErr.s1Rx*100),"%"));
if width(TTback) == 12
    disp(strcat("backErr.s2Rx = ",num2str(backErr.s2Rx*100),"%"));
end

% Plot the interpolated data
TTinterpF = fixerrorns(TTfront);
TTinterpB = fixerrors(TTback);

plotAurora(TTinterpF,strcat("FRONT ", name, ": Channel measurements - Interpolated"));
plotAurora(TTinterpB,strcat("BACK ", name, ": Channel measurements - Interpolated"));

%% Clip Data Set A
% Using local minima to find the front start point - where the motion
% begins, after the initial position data.

% Peak('P') and Trough ('T'): Values ('V') and Locations ('L')
[PVbutterF,PLbutterF,TVbutterF,TLbutterF] = multiPT(TTbuttercF,strcat("Front ", name, ": Filtered - Channel peaks and troughs"));

% Plot one singular channel for prettiness
figure('Name',strcat("Front ", name, ": Filtered - s1Tz peaks and troughs"));
findpeaks(TTbuttercF.s1Tz)
title(strcat("Front ", name, " Filtered - s1Tz Peaks and Troughs"))
text(PLbutterF.s1Tz+12,PVbutterF.s1Tz+.5,num2str((1:numel(PVbutterF.s1Tz))'))
text(TLbutterF.s1Tz+12,gnegate(TVbutterF.s1Tz)-.5,num2str((1:numel(TVbutterF.s1Tz))'))
hold on;
scatter(TLbutterF.s1Tz,gnegate(TVbutterF.s1Tz)-1,'^',"filled")

%% Edit here to find the trough that corresponds to the motion starting?
% take average of the two main "large" motion directions - s1Tx and
% s1Tz, at the "person-identified" trough #, then clip and save as
% "a" version ("b" version is back end, IF there aren't a billion
% missed data points. Less than 5%? Will adjust as I review all the
% datasets.)

% Indicate which # (and peak or trough)
% Exp1-1: T2,3   Exp1-2: T2     Exp1-3: T3 
% Exp2-1: T6,5   Exp2-2: T2,5 
% Exp5-1: T2     Exp5-2: T2,3   Exp5-3: T2 
% Exp6-1: T5     Exp6-2: T5     Exp6-3: T4

switch name
    case {'exp1-1','exp5-2'}
        ptStart = 2;
        ptStart2 = 3;
    case {'exp1-2', 'exp5-1', 'exp5-3'}
        ptStart = 2;
        ptStart2 = 2;
    case {'exp1-3'}
        ptStart = 3;
        ptStart2 = 3;
    case {'exp2-trigger-1'}
        ptStart = 6;
        ptStart2 = 5;
    case {'exp2-trigger-2'}
        ptStart = 2;
        ptStart2 = 5;
    case {'exp6-trigger-1', 'exp6-trigger-2'}
        ptStart = 5;
        ptStart2 = 5;
    case {'exp6-trigger-3'}
        ptStart = 4;
        ptStart2 = 4;
    otherwise
end

% For trough clips:
clipLocation = round(((TLbutterF.s1Tx(ptStart)+TLbutterF.s1Tz(ptStart2))/2),0,tieBreaker="even");

TTinit = TTbuttercF(1:clipLocation-1,:);
%writetimetable(TTinit,strcat('Datasets/aurora-',name,'-init.csv'));
TTinitR = TTfront(1:clipLocation-1,:);
%writetimetable(TTinitR,strcat('Datasets/Raw/aurora-',name,'-raw-init.csv'));

TTdataA = TTbuttercF(clipLocation:end,:);
%writetimetable(TTdataA,strcat('Datasets/aurora-',name,'-a.csv'));
TTdataRA = TTfront(clipLocation:end,:);
%writetimetable(TTdataRA,strcat('Datasets/Raw/aurora-',name,'-raw-a.csv'));

figure("Name",strcat(name, " Dataset 'A'"))
plot(TTdataA.s1Tz)
title(strcat(name, " Dataset 'A'"))

%% Dataset 'A' Error Percent

AErr = checkError(TTdataRA);
disp(strcat("AErr.s1Rx = ",num2str(AErr.s1Rx*100),"%"));
if width(TTdataRA) == 12
    disp(strcat("AErr.s2Rx = ",num2str(AErr.s2Rx*100),"%"));
end
plotAurora(TTdataRA,strcat(name, " Dataset 'A': Raw Channel data"));

%% Clip Data Set B
% Using local minima to find the back end point - where the motion
% ends.

% Peak('P') and Trough ('T'): Values ('V') and Locations ('L')
[PVbutterB,PLbutterB,TVbutterB,TLbutterB] = multiPT(TTbuttercB,strcat("Back ", name, ": Filtered - Channel peaks and troughs"));

% Plot one singular channel for prettiness
figure('Name',strcat("Back ", name, ": Filtered - s1Tz peaks and troughs"));
findpeaks(TTbuttercB.s1Tz)
title(strcat("Back ", name, " Filtered - s1Tz Peaks and Troughs"))
text(PLbutterB.s1Tz+12,PVbutterB.s1Tz+.5,num2str((1:numel(PVbutterB.s1Tz))'))
text(TLbutterB.s1Tz+12,gnegate(TVbutterB.s1Tz)-.5,num2str((1:numel(TVbutterB.s1Tz))'))
hold on;
scatter(TLbutterB.s1Tz,gnegate(TVbutterB.s1Tz),'^',"filled")

%% Edit here to find the trough that corresponds to the motion ending?
% take average of the two main "large" motion directions - s1Tx and
% s1Tz, at the "person-identified" trough #, then clip and save as
% "b" version (considering 5% to be error threshold)

% Indicate which # (and peak or trough)
% Exp1-1: P3     Exp1-2: P3     Exp1-3: P3 
% Exp2-1: P3     Exp2-2: P4 
% Exp5-1: P3     Exp5-2: P3     Exp5-3: P3
% Exp6-1: P3     Exp6-2: P3     Exp6-3: P4,3

switch name
    case {'exp1-1', 'exp1-2', 'exp1-3', 'exp2-trigger-1', 'exp5-1', ...
            'exp5-2', 'exp5-3', 'exp6-trigger-1', 'exp6-trigger-2'}
        ptEnd = 3;
        ptEnd2 = 3;
    case {'exp2-trigger-2'}
        ptEnd = 4;
        ptEnd2 = 4;
    case {'exp6-trigger-3'}
        ptEnd = 4;
        ptEnd2 = 3;
    otherwise
end

% For trough clips:
clipLocation = round(((PLbutterB.s1Tx(ptEnd)+PLbutterB.s1Tz(ptEnd2))/2),0,tieBreaker="even");

TTdataB = TTbuttercB(1:clipLocation-1,:);
%writetimetable(TTdataB,strcat('Datasets/aurora-',name,'-b.csv'));
TTdataRB = TTback(1:clipLocation-1,:);
%writetimetable(TTdataRB,strcat('Datasets/Raw/aurora-',name,'-raw-b.csv'));

TTlag = TTbuttercB(clipLocation:end,:);
%writetimetable(TTlag,strcat('Datasets/aurora-',name,'-lag.csv'));
TTlagR = TTback(clipLocation:end,:);
%writetimetable(TTlagR,strcat('Datasets/Raw/aurora-',name,'-raw-lag.csv'));

figure("Name",strcat(name, " Dataset 'B'"))
plot(TTdataB.s1Tz)
title(strcat(name, " Dataset 'B'"))

%% Dataset 'B' Error Percent

BErr = checkError(TTdataRB);
disp(strcat("BErr.s1Rx = ",num2str(BErr.s1Rx*100),"%"));
if width(TTdataRB) == 12
    disp(strcat("BErr.s2Rx = ",num2str(BErr.s2Rx*100),"%"));
end
plotAurora(TTdataRB,strcat(name, " Dataset 'B': Raw Channel data"));

%% Save Error Percents in Results

%clearvars -except name AErr BErr frontErr backErr
load 'Results\error.mat'

cName = [name(1:4),name(end)];

errorAll.(cName) = struct('AErr', AErr, 'BErr', BErr, 'frontErr', frontErr, 'backErr', backErr);
error.(strcat(cName,'_A_s1Rx')) = AErr.s1Rx;
error.(strcat(cName,'_A_s2Rx')) = AErr.s2Rx;
error.(strcat(cName,'_B_s1Rx')) = BErr.s1Rx;
error.(strcat(cName,'_B_s2Rx')) = BErr.s2Rx;

%save Results\error.mat errorAll error -mat