%{

atiPreProcess.m
    Import ATI data and filter. Plot Raw and Filtered channels.
    Clip at the "motion start" - the through at the beginning of the motion
        in Channel Fy. This corresponds to the start of the motion in the
        Aurora data. 
    Plot the peaks and troughs, and the final font-clipped dataset.
    Save the data that was trimmed from the front as "initial" data.

Save clipped (and init) data as raw:
    (Datasets/Raw: force-[experiment]-raw[-init].csv)
and filtered:
    (Datasets: force-[experiment][-init].csv)

Nicole Devos for the WearME Lab, Western University
Version 1.2
July 20, 2023

%}

close all
clear all
clc

name = "exp6-trigger-3";

lowhz = 0.5;
atiFreq = 1000/16;
VariableNames = [{"Fx"},{"Mx"},{"Fy"},{"My"},{"Fz"},{"Mz"}];
%rrrrrrrrrrrrrr4eeeeeeeeeeeeeeeeeegf2
TTati = readtimetable(strcat('FEB 13 DATA/force-feb-13-',name,'.csv'),'SampleRate',atiFreq,'Range','A:F');
TTati.Properties.VariableNames = [VariableNames{[1;3;5;2;4;6]}];
TTati.Time = duration(TTati.Time,'Format','hh:mm:ss.SSSS');

% Filter the data using a Butterworth lowpass filter, cutoff 18 hz, human
% motion  is about 2 hz, but try 4 hz?:

TTatiF = filterTT(TTati,atiFreq,lowhz);

disp(strcat('ATI data for ''', name, ''' imported'))

figure('Name',strcat(name,': Force and Moment data - Filtered'));
for i=1:length(VariableNames)
    subplot(3,2,i);
    plot(TTati.Time,TTati.(VariableNames{i}))
    hold on;
    plot(TTatiF.Time,TTatiF.(VariableNames{i}))
    title(VariableNames{i})
end

%% Clip at motion start

% Using local minima/maxima to find the front start point

% Peak('P') and Trough ('T'): Values ('V') and Locations ('L')
[PVatiF,PLatiF,TVatiF,TLatiF] = multiPT(TTatiF,strcat("Peaks and Troughs in the ", name," ATI data"));

% Plot one singular channel for prettiness
figure('Name',strcat("Peaks and troughs in the ", name, " Fy data"));
findpeaks(TTatiF.Fy)
title(strcat("Peaks and troughs in the ", name, " Fy data"))
text(PLatiF.Fy+12,PVatiF.Fy+.5,num2str((1:numel(PVatiF.Fy))'))
text(TLatiF.Fy+12,gnegate(TVatiF.Fy)-.5,num2str((1:numel(TVatiF.Fy))'))
hold on;
scatter(TLatiF.Fy,gnegate(TVatiF.Fy),'^',"filled")

%% Clip using one of the troughs
% Should be the same as the second number in auroraPreProcess.m
% 1-1: 3

% Trough location for clipping
switch name
    case {'exp1-1'}
        ptStart = 3;
    case {'exp6-trigger-3'}
        ptStart = 4;
    otherwise
end

clipLocation = TLatiF.Fy(ptStart);
    
TTatiInit = TTatiF(1:clipLocation-1,:);
%writetimetable(TTatiInit,strcat('Datasets/force-',name,'-init.csv'));
TTatiRI = TTati(1:clipLocation-1,:);
%writetimetable(TTatiRI,strcat('Datasets/Raw/force-',name,'-raw-init.csv'));

TTatiData = TTatiF(clipLocation:end,:);
%writetimetable(TTatiData,strcat('Datasets/force-',name,'.csv'));
TTatiRD = TTati(clipLocation:end,:);
%writetimetable(TTatiRD,strcat('Datasets/Raw/force-',name,'-raw.csv'));

figure("Name",strcat("Fy data from ", name));
plot(TTatiData.Fy)
title(strcat(name, ": Fy Data"))