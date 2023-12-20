%{

LOADDATA load data from mat files
MAT files created using importData.m, 

Nicole Devos for the WearME Lab, Western University

Version 1.0
January 4, 2023

%}

%% Measured Landmarks

% Measurement data from SW models. Measures from the landmark to the
% object's own COM. Note different base frame orientations.

% These might not be needed anymore
modelSkullLandmarkToCOM = [0;0;0];
cageLandmarkToCOM = [0;0;0];
wheelLandmarkToCOM = [0;0;0];

% auroraTo[]Landmark: Measurement data from Aurora. Measures the location
% of the Landmarks wrt the Aurora. This is just translation data, no
% rotation (Tx/y/z, no R).
% auroraToBase: Location of the base/ATI measured with the Aurora.
% jt2: location of joint 2's rotation axis (when in neutral position, as
% measured by the auopra)
% Calculated in importData.m
load MAT-Files/feb-13-distance-data.mat

% These might not be needed anymore
auroraToCOM = auroraToSkullLandmark + modelSkullLandmarkToCOM;
BaseToCOM = -auroraToBase + auroraToCOM;

auroraToCageCOM = auroraToCageLandmark + cageLandmarkToCOM; 
BaseToCageCOM = -auroraToBase + auroraToCageCOM;

auroraToWheelCOM = auroraToWheelLandmark + wheelLandmarkToCOM; 
BaseToWheelCOM = -auroraToBase + auroraToWheelCOM;


%% Import experiment data

% Import the specific aurora/ati dataset ('name')
load(strcat('MAT-Files/',name,'.mat'))

% Plot the data to confirm (plotAurora(TTraw, figure name, TTfiltered))
plotAurora(TTinterp, strcat("Aurora comparison between the raw and filtered data in the ", name, " experiment set"), TTa);
plotATI(TTatiRaw, strcat("ATI comparison between the raw and filtered data in the ", name, " experiment set"), TTati);

%% Initial position data

[sensor1, sensor2, plotInitAvgs] = initAvg(TTinit);
plotAurora(TTinit, strcat("Initial data vs average of first 100 points in the ", name, " dataset"),plotInitAvgs);