% Load constants

% Ver 1.0 April 17, 2023

%% DATA FROM SW

% Mass Data, in grams:
mh = 2384.53;%/1000; % head, cap, and three screws only
% mass of motor: 70g
mg = 113.60;%/1000; % cage, mounting hardware, motor, wheel hardware, arms.
% Measured in Dr. Price's lab: 56.0638 g
% Note: pre-2023: SW gave approximate weight as 55.70 g (??)
mw = 55.47;%/1000; % wheel only

% Orientation:
% x horizontal: +ve -> Skull left
% y vertical: +ve -> Crown
% z in/out: +ve -> Face

% Rotation matrix for SW->Ch4
% -ve quarter rotation about z, then -ve quarter rotation about y
Rswb = Rz(-pi/2)*Ry(-pi/2);

% In mm:
sw_SkullLandmarkToCOM = [1.33;48.24;-106.04];
sw_S2ToCageCOM = [-164.00;-46.54;-0.08];
sw_S2ToWheelCOM = [-164.00;-28.96;-0.08];

% Convert to Base frame:
skullLandmarkToCOM = transpose(Rswb)*sw_SkullLandmarkToCOM;%./1000;
s2ToCageCOM = transpose(Rswb)*sw_S2ToCageCOM;%./1000;
s2ToWheelCOM = transpose(Rswb)*sw_S2ToWheelCOM;%./1000;

% Axes of inertia, at COM, grams*mm^2:
% Principle axes
ih = [0; 0; 1];
jh = [0; -1; 0];
kh = [1; 0; 0];
uh = [ih, jh, kh];
iw = [0; 1; 0];
jw = [0; 0; 1];
kw = [1; 0; 0];
uw = [iw, jw, kw];

% Principle inertias in SW frame:
sw_IH = [6146383.34; 8123580.47; 9719428.66];%./1000;
sw_IH = abs(uh*sw_IH);
sw_IG = [68809.34; 87013.14; 142392.58];%./1000;
sw_IW = [13835.29; 19582.42; 19582.42];%./1000;
sw_IW = abs(uw*sw_IW);

% Convert to Base Frame:
IH = transpose(Rswb)*sw_IH;
IG = transpose(Rswb)*sw_IG;
IW = transpose(Rswb)*sw_IW;

% set Ixx, Iyy, and Izz for ease of reading code:
IHxx = IH(1);
IHyy = IH(2);
IHzz = IH(3);
IGxx = IG(1);
IGyy = IG(2);
IGzz = IG(3);
IWxx = IW(1);
IWyy = IW(2);
IWzz = IW(3);

save MAT-Files/gyro-consts.mat skullLandmarkToCOM s2ToCageCOM s2ToWheelCOM mh mg mw IHxx IHyy IHzz IGxx IGyy IGzz IWxx IWyy IWzz
%save MAT-Files/gyro-consts_m.mat skullLandmarkToCOM s2ToCageCOM s2ToWheelCOM mh mg mw IHxx IHyy IHzz IGxx IGyy IGzz IWxx IWyy IWzz

%% NOTES
% mneck = 0;
% mhead = 0;
% mcap = 0;
% mscrew = 0;
% mheadarm = 0;
% mcagebase = 0;
% mpanel = 0;
% mmotor = 0;
% mcage = 0;
% mcagearm = 0;
% mwheel = 0;
% mnut = 0;
% 
% mh = mneck+mhead+mcap+4*mscrew+mheadarm+mcagebase+4*mpanel;
% mg = mmotor+mcage+2*mcagearm;
% mw = mwheel+2*mnut;
% 90[;p'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''