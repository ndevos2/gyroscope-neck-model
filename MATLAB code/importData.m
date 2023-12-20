%{

IMPORTDATA Import, pre-calculate (where needed), and save variables as mat
files. Run once.

Nicole Devos for the WearMe Lab, Western University

Version 2.1
July 20, 2023

%}

%% 

close all
clear variables

%% Import Landmark location wrt Aurora
% These might not be needed anymore

auroraToSkullLandmark = landmarkAvg('FEB 13 CLEAN DATA/aurora-clean-feb-13-skull-landmark.csv')./1000;
auroraToCageLandmark = landmarkAvg('FEB 13 CLEAN DATA/aurora-clean-feb-13-cage-landmark.csv')./1000;
auroraToWheelLandmark = landmarkAvg('FEB 13 CLEAN DATA/aurora-clean-feb-13-wheel-landmark.csv')./1000;

%% Import Base location wrt Aurora

auroraToBase = [positionAvg('FEB 13 CLEAN DATA/aurora-clean-feb-13-base-x-z.csv','E:E');positionAvg('FEB 13 CLEAN DATA/aurora-clean-feb-13-base-y.csv','F:F');positionAvg('FEB 13 CLEAN DATA/aurora-clean-feb-13-base-x-z.csv','G:G')]./1000;

%% Import Joint 2 location

% Note: height = x-axis: 1
jt2 = landmarkAvg('FEB 13 CLEAN DATA/aurora-clean-feb-13-jt2-x')./1000;

%% Save above variables

save MAT-Files\feb-13-distance-data_m.mat

clear all

%% Load and save raw, filtered Aurora and ATI data to one MAT file.

% Load all file names; RA: Aurora, raw set A. RI: Aurora, raw init. A:
% Aurora, set A. I: Aurora, init. FR: Force, raw. FRI: Force, raw init. F:
% Force. FI: Force, init. B: Aurora, set B.

listRA = dir('Datasets\Raw\aurora-exp*-raw-a.csv');
listRB = dir('Datasets\Raw\aurora-exp*-raw-b.csv');
listRI = dir('Datasets\Raw\aurora-exp*-raw-init.csv');
listA = dir('Datasets\aurora-exp*-a.csv');
listB = dir('Datasets\aurora-exp*-b.csv');
listI = dir('Datasets\Init\aurora-exp*-init.csv');
listFR = dir('Datasets\Raw\force-exp*-raw.csv');
listFRI = dir('Datasets\Raw\aurora-exp*-raw-init.csv');
listF = dir('Datasets\force-exp*.csv');
listFI = dir('Datasets\Init\force-exp*-init.csv');

% loop through the files, importing the TTs and saving as workspace
% variables.
for i = 1:length(listRA)
    % Raw Aurora
    TTrawA = readtimetable(strcat('Datasets/Raw/',listRA(i).name));
    TTrawB = readtimetable(strcat('Datasets/Raw/',listRB(i).name));
    TTinterpA = fixerrors(TTrawA);
    TTinterpB = fixerrors(TTrawB);
    TTRinit = readtimetable(strcat('Datasets/Raw/',listRI(i).name));
    % Filtered Aurora
    TTa = readtimetable(strcat('Datasets/',listA(i).name));
    TTb = readtimetable(strcat('Datasets/',listB(i).name));
    TTinit = readtimetable(strcat('Datasets/Init/',listI(i).name));
    
    % Raw ATI (Force/Torque)
    TTatiRaw = readtimetable(strcat('Datasets/Raw/',listFR(i).name));
    TTatiRawinit = readtimetable(strcat('Datasets/Raw/',listFRI(i).name));
    % Filtered ATI (Force/Torque)
    TTati = readtimetable(strcat('Datasets/',listF(i).name));
    TTatiInit = readtimetable(strcat('Datasets/Init/',listFI(i).name));

    varSaveName = listF(i).name(7:end-4);
    save(strcat('MAT-Files/',varSaveName),'TTrawA','TTrawB','TTinterpA','TTinterpB','TTRinit','TTa','TTb','TTinit','TTatiRaw','TTatiRawinit','TTati','TTatiInit')
    clear TTrawA TTrawB TTinterpA TTinterpB TTRinit TTa TTb TTinit TTatiRaw TTatiRawinit TTati TTatiInit
end