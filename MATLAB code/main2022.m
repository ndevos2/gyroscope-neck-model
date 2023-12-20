%{

"DESCRIPTION"
Note:
sensor1 MUST be positioned on the head and sensor2 MUST be positioned
on the cage. ATI sensor must be zeroed before collecting data.

Reminder:
for this analysis: Left/Right is X, Front/Back (Anterior/Posterior)
is Y, and Up/Down (Superior/Inferior) is Z

Nicole Devos for the WearME Lab, Western University

Version 4.2
July 20, 2023
6yh66
%}

%% Initialization

close all
clearvars

name = 'exp6-trigger-3';
set = 'B';
% might not use these sets: exp1-1, 1-2, 1-3, 2-trigger-2
% I think this one failed the inclusion criteria: 2-trigger-1
% better: 5-1, 5-2, 5-3, 6-trigger-1, 6-trigger-2, 6-trigger-3

% sampling frequencies, in Hz
aurFreq = 40;
atiFreq = 1000/16;

% Device: motor speed (voltage * motor "kV")
voltage = 12;
kV = 950;

%% Data Import

%load data script 
load MAT-Files/feb-13-distance-data.mat
% Import the specific aurora/ati dataset ('name')
load(strcat('MAT-Files/',name,'.mat'))

%% Equation constants

% needs to be  in mm/s^2
g = 9.81*1000; % m/s^2 -> convert to mm/s^2 if needed

% Omega: rad*/s
% RPM = voltage*kV;
% rotation/min * 2*pi radians/rotation * 1/60 min/sec
Omega = voltage*kV*(2*pi)/60; 

% Created with loadConstants.m
load MAT-Files/gyro-consts.mat

% masses - compare real to predicted in model? Then can comment on the
% accuracy/inaccuracy of the predicted model mass moments of inertia.

% Calculate the distance from Base to COMs, system frame (N)
Rab = Rz(-pi/2)*Rx(pi/2);

dbH = ( transpose(Rab)*auroraToSkullLandmark + skullLandmarkToCOM ) - transpose(Rab)*auroraToBase;
dbG = ( transpose(Rab)*auroraToCageLandmark + s2ToCageCOM ) - transpose(Rab)*auroraToBase;
dbW = ( transpose(Rab)*auroraToWheelLandmark + s2ToWheelCOM ) - transpose(Rab)*auroraToBase;

xh = dbH(1);
yh = dbH(2);
zh = dbH(3);

xg = dbG(1);
yg = dbG(2);
zg = dbG(3);

xw = dbW(1);
yw = dbW(2);
zw = dbW(3);

%% Aurora Analysis (Chapter 6)

% Rotation matrix for aurora (A) to base frame (0)
RA0 = Rz(pi)*Rx(pi/2);
dA0 = auroraToBase; % position of Base Frame, measured with Aurora

% inverse (to convert Aurora data to Base frame) - confirmed with animation
H0A = [transpose(RA0), -transpose(RA0)*dA0; 0 0 0 1];

% das1 is starting position of sensor. d0s1 is position of sensor1 wrt Base
das1 = mean(TTinit{1:50,4:6});%./1000;
das2 = mean(TTinit{1:50,10:12});%./1000;

% Testing deviation in initial position(s)
%rmse50a = sqrt(mean((TTinit{1:50,4:6} - das1).^2));
%das2 = mean(TTinit{1:50,10:12});
%rmse50b = sqrt(mean((TTinit{1:50,10:12} - das2).^2));

% [l1 l3 l2] = d0s1 = Ra0*(das1-dab);
d0s1 = transpose(RA0)*(das1'-auroraToBase);
l1 = abs(d0s1(1));
l2 = abs(d0s1(3));
l3 = abs(d0s1(2));

% l4 = distance from base to joint 2 along y0 axis
% [l4 yc0 l5] = d0c = Ra0*((dacl+dclcx)-dab);
d0c = transpose(RA0)*(jt2-auroraToBase);
l4 = abs(d0c(2));

% [x0s2 l7 l5] = d0s2 = Ra0*(das2-dab);
d0s2 = transpose(RA0)*(das2'-auroraToBase);
l5 = abs(d0s2(3));
l6 = abs(d0s2(1))-l4;
l7 = abs(d0s2(2));

% empty matrices - TTa/b is the Aurora data timetable
switch set
    case 'A'
        TTset = TTa;
    case 'B'
        TTset = TTb;
        otherwise%qawssssssssssss3444
end
pos_s1 = zeros(height(TTset),3);
pos_s2 = zeros(height(TTset),3);
convertedAngles = zeros(height(TTset),2);
calculatedAngles = zeros(height(TTset),2);

% calculate sensor position and angle for each data point
for i=1:height(TTset)
    Ras1 = eulerConv(TTset{i,1:3}); % convert Aurora data (euler angles) to Rotation Matrix
    das1 = TTset{i,4:6}';%./1000; % isolate position data from the TT
    Tas1 = [Ras1,das1;0 0 0 1]; % T matrix of s1 position
    T02_s1 = H0A*Tas1; % T02 for s1: includes end effector frame
    
    % save position xyz in a matrix
    pos_s1(i,:) = T02_s1(1:3,4)';
    
    % Repeat above for second sensor data:
    Ras2 = eulerConv(TTset{i,7:9});
    das2 = TTset{i,10:12}';%./1000;
    Tas2 = [Ras2 das2; 0 0 0 1];
    T02_s2 = H0A*Tas2;

    % save position xyz in a matrix
    pos_s2(i,:) = T02_s2(1:3,4)';
    
    % for ease of reading:
    x02_1 = T02_s1(1,4);
    y02_1 = T02_s1(2,4);
    x02_2 = T02_s2(1,4);
    y02_2 = T02_s2(2,4);
    z02_2 = T02_s2(3,4);

    % calculate c1, s1; c2, s2
    c1 = (x02_1*l1+y02_1*l3)/(l1^2+l3^2);
    s1 = (y02_1 - l3*c1)/l1;
%    s2 = -(z02_2+l5)/l6;
%    c2 = (l4^2 + 3*l5^2 + 2*l5*z02_2 + l6^2 + l7^2 - x02_2^2 - y02_2^2 - z02_2^2)/(2*l4*l6);
    s2 = (-z02_2+l5)/l6;
    %c2 = (x02_2^2 + y02_2^2 + z02_2^2 - (l4^2 + l5^2 + l6^2 + l7^2) + 2*l5*l6*s2)/(2*l4*l6);
    c2 = (x02_2^2 + y02_2^2 + z02_2^2 - (l4^2 + l5^2 + l6^2 + l7^2) + 2*l5^2 + 2*l5*z02_2)/(2*l4*l6);
    %using x/y arent quite the same
    %c2 = (-y02_2+l4*s1+l7*c1)/(l6*s1);
    
    convertedAngles(i,:) = [atan2d(s1,c1),atan2d(s2,c2)];
    calculatedAngles(i,:) = [atan2(s1,c1),atan2(s2,c2)];
end

% Calculated angle is the opposite direction of desired angle
% Only theta, head angle?
convertedAngles(:,1) = gnegate(convertedAngles(:,1));
calculatedAngles(:,1) = gnegate(calculatedAngles(:,1));

figure('Name',"phi in degrees");
plot(convertedAngles(:,1))
title(strcat("phi (head angle) in degrees for ",name),"Presented in the Aurora conversion frame")

figure('Name',"alpha in degrees");
plot(convertedAngles(:,2))
title(strcat("alpha (cage angle) in degrees for ",name),"Presented in the Aurora conversion frame")
%% Wrap check for Cage angles
%{
for i=1:length(convertedAngles)
   if convertedAngles(i,2) < 0
       convertedAngles(i,2) = convertedAngles(i,2) + 360;
   end
end

figure;
plot(convertedAngles(:,2))

figure;
sgtitle("calculated head angle, cage angle");
subplot(2,1,1);
plot(convertedAngles(:,1));
title("head angle")
hold on;
subplot(2,1,2);
plot(convertedAngles(:,2));
title("cage angle")
hold off;
%}

%% Position, Velocity, Acceleration

% Velocity
calculatedSpeed = zeros(length(calculatedAngles)-1,2);

calculatedSpeed(:,1) = diff(calculatedAngles(:,1))/(1/aurFreq);
calculatedSpeed(:,2) = diff(calculatedAngles(:,2))/(1/aurFreq);

% Acceleration
calculatedAccel = zeros(length(calculatedAngles)-2,2);

calculatedAccel(:,1) = diff(calculatedSpeed(:,1))/(1/aurFreq);
calculatedAccel(:,2) = diff(calculatedSpeed(:,2))/(1/aurFreq);

% Save to Timetable

TTmotion = array2timetable([convertedAngles(3:end,1),calculatedAngles(3:end,1),calculatedSpeed(2:end,1),calculatedAccel(:,1),convertedAngles(3:end,2),calculatedAngles(3:end,2),calculatedSpeed(2:end,2),calculatedAccel(:,2)],'SampleRate',aurFreq);
TTmotion.Properties.VariableNames = ["deg1","pos1","vel1","accel1","deg2","pos2","vel2","accel2"];
save(strcat("Mat-Files\",name,"-TTmotion"),'TTmotion');

% Plot
figure('Name',strcat(name,': head and cage motion data'));
sgtitle("Head and cage motion data")
subplot(3,2,1)
plot(TTmotion.Time, TTmotion.pos1);
title('Head position')
ylabel('rad')
hold on;

subplot(3,2,2)
plot(TTmotion.Time, TTmotion.pos2);
title('Cage position')
ylabel('rad')

subplot(3,2,3)
plot(TTmotion.Time, TTmotion.vel1);
title('Head velocity')
ylabel('rad/s')

subplot(3,2,4)
plot(TTmotion.Time, TTmotion.vel2);
title('Cage velocity')
ylabel('rad/s')

subplot(3,2,5)
plot(TTmotion.Time, TTmotion.accel1);
title('Head acceleration')
ylabel('rad/s^2')

subplot(3,2,6)
plot(TTmotion.Time, TTmotion.accel2);
title('Cage acceleration')
ylabel('rad/s^2')
hold off;

%% Synchronize/interpolation for motion and force data

% Re-time TTati to start at 0sec
TTati = table2timetable(removevars(timetable2table(TTati),"Time"),'SampleRate',atiFreq);
TTa = table2timetable(removevars(timetable2table(TTa),"Time"),'SampleRate',aurFreq);

switch set
    case 'A'
        TTati = TTati(TTati.Time <= TTa.Time(end-2), :); % First half ATI
        TTset = TTa;
    case 'B'
        TTb = table2timetable(removevars(timetable2table(TTb),"Time"),'SampleRate',aurFreq);
        TTati = TTati((TTati.Time > TTa.Time(end-2)), :); % Second half ATI
        TTati = table2timetable(removevars(timetable2table(TTati),"Time"),'SampleRate',atiFreq);
        TTati = TTati((TTati.Time <= TTb.Time(end-2)), :); % Clip trailing data
        TTset = TTb;
    otherwise
        TTati = 0;
end


TTsynch = synchronize(TTmotion,TTati,'union','linear');

% Head motion: phi
phi = TTsynch.pos1;
dphi = TTsynch.vel1;
ddphi = TTsynch.accel1;

% Cage motion: alpha
alpha = TTsynch.pos2;
dalpha = TTsynch.vel2;
ddalpha = TTsynch.accel2;

%%%%%%%%%%% NEED TO MOVE THE ABOVE INTO CH4 COORDINATES %%%%%%%%%%%%

%% Theoretical Model (Chapter 4)

k1 = mh*(zh^2 + yh^2) + (mg + mw)*(zg^2 + yg^2) + IHxx;
k2 = IGxx + IWxx;
k3 = IGzz + IWzz;
k4 = g*(mh*zh + (mg + mw)*zg);
k5 = g*(mh*yh + (mg + mw)*yg);
k6 = IGyy + IWxx;

M1 = (k1 + k2).*ddphi + (k3 - k2).*sin(alpha).^2.*ddphi + (k3 - k2).*sin(2.*alpha).*dalpha.*dphi + IWxx.*Omega.*cos(alpha).*dalpha - k4.*sin(phi) + k5.*cos(phi);
M2 = k6.*ddalpha - (k3 - k2).*sin(alpha).*cos(alpha).*dphi.^2 - IWzz.*Omega.*cos(alpha).*dphi;

% Convert to N and Nm
M1 = M1./1000./1000000; % 1kg/1000g and 1m^2/1,000,000 mm^2
M2 = M2./1000./1000000; % 1kg/1000g and 1m^2/1,000,000 mm^2

% Re-add bias to the sensor data
bias = mh*g*yh/1000/1000000+(mg+mw)*g*yg/1000/1000000;
Mbias = gnegate(TTsynch.Mx)+bias;

%% Comparison, Results

% No bias
figure('Name',"Comparing M1 calculated and measured");
plot(TTsynch.Time,M1)
hold on;
plot(TTsynch.Time,gnegate(TTsynch.Mx),'r')
title(strcat("Comparing Calculated and Measured Moment for ", name))
legend('theoretical torque (Nm)','measured torque (Nm)','Location','northwest')
% NW for exp1-1A
% N for exp1-1B ->
% S for exp21B

% Add in bias:
figure('Name',"Comparing M1 and actual (test)");
plot(TTsynch.Time,M1)
hold on;
plot(TTsynch.Time,Mbias,'r')
title(strcat("Comparing Calculated and (test) Actual Moment for ", name))
legend('theoretical torque (Nm)','measured torque (Nm)','Location','southwest')

% Adjusting for the human error in synchronizing the datasets:
% step 1: make new tables
TTactual = TTsynch(:,"Mx");
TTactual.Mx = Mbias;
TTactual = renamevars(TTactual,"Mx","Mbias");

TTcalc = TTsynch(:,"accel1");
TTcalc.accel1 = M1;
TTcalc = renamevars(TTcalc,"accel1","M1");

% step 2: find local minima & maxima
actualMinBool = islocalmin(TTactual.Mbias);
actualMaxBool = islocalmax(TTactual.Mbias);
calcMinBool = islocalmin(TTcalc.M1);
calcMaxBool = islocalmax(TTcalc.M1);

actualExtremaBool = actualMinBool|actualMaxBool;
calcExtremaBool = calcMinBool|calcMaxBool;

% Plot
figure('name',"Using local minima & maxima to adjust time offset")
plot(TTactual.Time,TTcalc.M1,TTactual.Time,TTactual.Mbias)
title("Using local minima & maxima to adjust for time offset",name)
hold on;
plot(TTactual.Time(calcExtremaBool),TTcalc.M1(calcExtremaBool),'.',TTactual.Time(actualExtremaBool),TTactual.Mbias(actualExtremaBool),'*')
legend('theoretical torque (Nm)','measured torque (Nm)','theoretical extrema','measured extrema','location','northwest');
%SE for exp11A
%S for exp11B

%%
% Remove leading or lagging extrema if necessary:
switch set
    case 'A'
    switch name
        case {'exp5-1','exp6-trigger-3'}
        case 'exp1-3'
            ind = find(calcExtremaBool);
            calcExtremaBool(ind(1)) = false;
        case 'exp5-2'
            ind = find(calcExtremaBool);
            calcExtremaBool(ind(1)) = false;
            ind = find(actualExtremaBool);
            actualExtremaBool(ind(1)) = false;
        otherwise
            ind = find(actualExtremaBool);
            actualExtremaBool(ind(1)) = false;
    end
    case 'B'
        switch name
            case {'exp1-1','exp5-1','exp6-trigger-1'}
            case 'exp6-trigger-3'
                %ind = find(actualExtremaBool);
                %actualExtremaBool(ind(end-1:end)) = false;
            otherwise
                ind = find(calcExtremaBool);
                calcExtremaBool(ind(end)) = false;
        end
end

% step 3: find time difference
actualMin = find(actualExtremaBool);
calcMin = find(calcExtremaBool);

offsetSet = TTactual.Time(actualMin)-TTcalc.Time(calcMin);

% UPDATE FOR EACH ONE
switch set
    case 'A'
        switch name
            case 'exp1-1'
                offset = seconds(0.088); % [0.088 0.060 0.062 0.070 0.071] [0.063, 0.085, 0.162, 0.070, 0.171]
            case 'exp1-2'
                offset = seconds(0.222); % [0.222 0.193 0.168 0.217 0.181] [0.247, 0.193, 0.218, 0.192, 0.206]
            case 'exp1-3'
                offset = seconds(0.026); %[0.026, -0.02, -0.002, -0.005, -0.014] [0.051, -0.02, -0.002, -0.005, -0.014]
            case 'exp2-trigger-1'
                offset = seconds(-0.005); % [-0.005 0.007 0.010 -0.013 0.011] [-0.030, 0.032, 0.035, -0.013, 0.111]
            case 'exp2-trigger-2'
                offset = seconds(0.061); % [0.061 0.030 0.046 0.035 0.040 0.036 0.038 0.057 0.051 0.047 0.031 0.056 0.022 0.036 0.024] [0.136, 0.055, -0.054, 0.06, -0.01, 0.061, 0.188, 0.057, 0.051, 0.047, -0.119, 0.031, 0.047, 0.136, 0.124]
            case 'exp5-1'
                offset = seconds(0.063); % [0.063 0.072 0.064 0.065 0.065] [0.088, 0.047, 0.114, 0.065, 0.065]
            case 'exp5-2'
                offset = seconds(0.025); % [0.025 0.004 0.029 0.015 0.005] [0.050, 0.004, 0.004, 0.015, 0.030]
            case 'exp5-3'
                offset = seconds(0.036); % [0.036 0.024 0.039 0.030 0.039] [0.061, 0.024, 0.039, 0.030, 0.039]
            case 'exp6-trigger-1'
                offset = seconds(0.026); % [0.026 0.027 0.015 0.031 0.028] [0.026, 0.002, 0.040, 0.031, -0.022]
            case 'exp6-trigger-2'
                offset = seconds(-0.004); % [-0.004 -0.006 -0.001 0.013 0.004] [0.021, 0.019, 0.024, 0.013, -0.021]
            case 'exp6-trigger-3'
                offset = seconds(0.005); % [0.005 -0.008 -0.019 -0.027 -0.009] [0.005, -0.008, 0.006, -0.052, -0.034]
            otherwise
                offset = seconds(100);
        end
    case 'B'
        switch name
            case 'exp1-1'
                offset = seconds(0.095); % [0.095 0.056 0.060 0.080 0.050] [0.147 0.110 0.115 0.106 0.104] [0.097 0.135 0.215 0.106 0.129]
            case 'exp1-2'
                offset = seconds(0.262); % [0.262 0.241 0.245 0.233 0.244] [ 0.262 0.216 0.27 0.233 0.244 ]
            case 'exp1-3'
                offset = seconds(0.082); % [0.082 0.044 0.042 0.055 0.041] [ 0.082 0.069 0.0442 0.055 0.41 ]
            case 'exp2-trigger-1'
                offset = seconds(0.145); % [0.145 0.083 0.062 0.073 0.057] [ 0.145 0.058 0.112 0.073 0.057 ]
            case 'exp2-trigger-2'
                offset = seconds(0.112); % [0.112 0.098 0.086 0.093 0.090 0.103 0.094] [ 0.087 0.073 0.111 0.118 0.090 0.103 0.069 ]
            case 'exp5-1'
                offset = seconds(0.132); % [0.132 0.114 0.116 0.118 0.127 0.010] [ 0.157 0.114 0.116 0.118 0.102 ]
            case 'exp5-2'
                offset = seconds(0.075); % [0.075 0.069 0.071 0.073 0.063] [ 0.125 0.069 0.046 0.073 0.113 ]
            case 'exp5-3'
                offset = seconds(0.101); % [0.101 0.094 0.080 0.084 0.085] [ 0.101 0.119 0.080 0.0840 0.060 ]
            case 'exp6-trigger-1'
                offset = seconds(0.097); % [0.097 0.099 0.092 0.101 0.085] [ 0.097 0.099 0.092 0.101 0.085 ]
            case 'exp6-trigger-2'
                offset = seconds(0.085); % [0.085 0.062 0.048 0.078 0.072] [ 0.085 0.062 0.73 0.003 0.072 ]
            case 'exp6-trigger-3'
                offset = seconds(0.095); % [0.095 0.056 0.060 0.080 0.050] [ 0.095 0.056 0.060 0.080 0.075 -0.044 -0.162 ]
            otherwise
                offset = seconds(100);
        end
end

% step 4: add difference
TTcalc.Time = TTcalc.Time+offset;

% Step 5: Re-plot
% Add in bias:
figure('Name',"Comparing M1 and actual, time adjusted");
plot(TTcalc.Time,TTcalc.M1)
hold on;
plot(TTactual.Time,TTactual.Mbias,'r')
title(strcat("Comparing Time-Adjusted Calculated and Actual Moment for ", name))
legend('theoretical torque (Nm)','measured torque (Nm)','Location','southwest')

% Step 6: Create new synchronized dataset
TTfixed = synchronize(TTactual,TTcalc,'union','linear');

%%  Analysis

% Root mean square error b/n forecast/predicted and actual

% Predicted
F = TTfixed.M1;
% Actual
A = TTfixed.Mbias;

% RMSD and then normalize
E = abs(A-F);
RMSE = sqrt(mean(E.^2));
RMSEn = sqrt(mean(E.^2))/(abs(max(TTfixed.Mbias)-min(TTfixed.Mbias)));
[maxE,maxEindex] = max(E);

RMSEunbiased = sqrt(mean(abs(gnegate(TTsynch.Mx)-M1).^2));
RMSEbiased = sqrt(mean(abs(Mbias-M1).^2));

figure('Name',"some kinda error test");
plot(TTfixed.Time,TTfixed.M1,TTfixed.Time,TTfixed.Mbias)
hold on;
title(strcat(name,": Actual vs Theoretical"), "Re-Biased and Time-Adjusted with maximum error shown")
plot([TTfixed.Time(maxEindex);TTfixed.Time(maxEindex)],[TTfixed.M1(maxEindex);TTfixed.Mbias(maxEindex)],'k-.')
legend('theoretical torque','actual torque','maximum error','Location','northwest')
% exp 21A NE 

%% Save to results folder

% load 'Results\results_fixed.mat'
% 
% cName = [name(1:4),name(end),set];
% 
% resultsAll.(cName) = struct('RMSE', RMSE, 'RMSEn', RMSEn, 'RMSEbiased', RMSEbiased, 'RMSEunbiased', RMSEunbiased, 'maxErr', maxE, 'maxActual', max(TTfixed.Mbias), 'minActual', min(TTfixed.Mbias));
% results.(strcat(cName,'_RMSEn')) = RMSEn;
% results.(strcat(cName,'_EMax')) = maxE;
% 
% save Results\results_fixed.mat resultsAll results -mat

TTunbiased = TTsynch(:,"Mx");

% Predicted
F = TTcalc.M1;
% Actual
A = TTunbiased.Mx;

% RMSD and then normalize
E = abs(A-F);
oriRMSEn = sqrt(mean(E.^2))/(abs(max(TTunbiased.Mx)-min(TTunbiased.Mx)));

load 'Results\results_fixed_normalized.mat'

cName = [name(1:4),name(end),set];

resultsNorm.(strcat(cName,'_RMSEn')) = RMSEn;
resultsNorm.(strcat(cName,'_oriRMSEn')) = oriRMSEn;

save Results\results_fixed_normalized.mat resultsNorm -mat

%% Just the non-interpolated ones
% % still need to interpolate to get the filtered data though...
% 
% switch set
%     case 'A'
%         TTrawS = TTrawA;
%     case 'B'
%         TTrawS = TTrawB;
%     otherwise
% end
% 
% % Re-time TTraw to start at 0:
% TTrawS = table2timetable(removevars(timetable2table(TTrawS),"Time"),'SampleRate',aurFreq);
% 
% TTset2 = TTset;
% TTrawS2 = TTrawS;
% sinput = TTrawS.Variables;
% soutput = TTset2.Variables;
% 
% % Replace error pts from filtered data:
% soutput(sinput<-1*10^10) = missing;
% 
% % Replace error pts from raw data:
% sinput(sinput<-1*10^10) = missing;
% 
% TTset2{:,:} = soutput;
% TTrawS2{:,:} = sinput;
% 
% plotAurora(TTset2,strcat(name, "-", set, ": Removing the error data points and comparing the Raw and Filtered Data"),TTrawS2);

% %% Re-Running Data
% 
% % empty matrices - TTa/b is the Aurora data timetable
% pos_s12 = zeros(height(TTset2),3);
% pos_s22 = zeros(height(TTset2),3);
% convertedAngles2 = zeros(height(TTset2),2);
% calculatedAngles2 = zeros(height(TTset2),2);
% 
% % calculate sensor position and angle for each data point
% for i=1:height(TTset2)
%     Ras1 = eulerConv(TTset2{i,1:3}); % convert Aurora data (euler angles) to Rotation Matrix
%     das1 = TTset2{i,4:6}';%./1000; % isolate position data from the TT
%     Tas1 = [Ras1,das1;0 0 0 1]; % T matrix of s1 position
%     T02_s1 = H0A*Tas1; % T02 for s1: includes end effector frame
%     
%     % save position xyz in a matrix
%     pos_s1(i,:) = T02_s1(1:3,4)';
%     
%     % Repeat above for second sensor data:
%     Ras2 = eulerConv(TTset2{i,7:9});
%     das2 = TTset2{i,10:12}';%./1000;
%     Tas2 = [Ras2 das2; 0 0 0 1];
%     T02_s2 = H0A*Tas2;
% 
%     % save position xyz in a matrix
%     pos_s2(i,:) = T02_s2(1:3,4)';
%     
%     % for ease of reading:
%     x02_1 = T02_s1(1,4);
%     y02_1 = T02_s1(2,4);
%     x02_2 = T02_s2(1,4);
%     y02_2 = T02_s2(2,4);
%     z02_2 = T02_s2(3,4);
% 
%     % calculate c1, s1; c2, s2
%     c1 = (x02_1*l1+y02_1*l3)/(l1^2+l3^2);
%     s1 = (y02_1 - l3*c1)/l1;
%     s2 = -(z02_2+l5)/l6;
%     c2 = (l4^2 + 3*l5^2 + 2*l5*z02_2 + l6^2 + l7^2 - x02_2^2 - y02_2^2 - z02_2^2)/(2*l4*l6); 
%     %using x/y arent quite the same
%     %c2 = (-y02_2+l4*s1+l7*c1)/(l6*s1);
%     
%     convertedAngles2(i,:) = [atan2d(s1,c1),atan2d(s2,c2)];
%     calculatedAngles2(i,:) = [atan2(s1,c1),atan2(s2,c2)];
% end
% 
% % Calculated angle is the opposite direction of desired angle
% convertedAngles2 = gnegate(convertedAngles2);
% calculatedAngles2 = gnegate(calculatedAngles2);
% 
% figure('Name',"phi in degrees");
% plot(convertedAngles(:,1))
% title(strcat("phi (head angle) in degrees for ",name),"Presented in the Aurora conversion frame")
% 
% figure('Name',"alpha in degrees");
% plot(convertedAngles2(:,2))
% title(strcat("alpha (cage angle) in degrees for ",name),"Presented in the Aurora conversion frame")
% %% Position, Velocity, Acceleration
% 
% % Velocity
% calculatedSpeed2 = zeros(length(calculatedAngles2)-1,2);
% 
% calculatedSpeed2(:,1) = diff(calculatedAngles2(:,1))/(1/aurFreq);
% calculatedSpeed2(:,2) = diff(calculatedAngles2(:,2))/(1/aurFreq);
% 
% % Acceleration
% calculatedAccel2 = zeros(length(calculatedAngles2)-2,2);
% 
% calculatedAccel2(:,1) = diff(calculatedSpeed2(:,1))/(1/aurFreq);
% calculatedAccel2(:,2) = diff(calculatedSpeed2(:,2))/(1/aurFreq);
% 
% % Save to Timetable
% 
% TTmotion2 = array2timetable([convertedAngles2(3:end,1),gnegate(calculatedAngles2(3:end,1)),calculatedSpeed2(2:end,1),calculatedAccel2(:,1),convertedAngles2(3:end,2),gnegate(calculatedAngles2(3:end,2)),calculatedSpeed2(2:end,2),calculatedAccel2(:,2)],'SampleRate',aurFreq);
% TTmotion2.Properties.VariableNames = ["deg1","pos1","vel1","accel1","deg2","pos2","vel2","accel2"];
% save(strcat("Mat-Files\",name,"-TTmotion2"),'TTmotion2');
% 
% % Plot
% figure('Name',strcat(name,': head and cage motion data'));
% sgtitle("Head and cage motion data")
% subplot(3,2,1)
% plot(TTmotion2.Time, TTmotion2.pos1);
% title('Head position')
% ylabel('rad')
% hold on;
% 
% subplot(3,2,2)
% plot(TTmotion2.Time, TTmotion2.pos2);
% title('Cage position')
% ylabel('rad')
% 
% subplot(3,2,3)
% plot(TTmotion2.Time, TTmotion2.vel1);
% title('Head velocity')
% ylabel('rad/s')
% 
% subplot(3,2,4)
% plot(TTmotion2.Time, TTmotion2.vel2);
% title('Cage velocity')
% ylabel('rad/s')
% 
% subplot(3,2,5)
% plot(TTmotion2.Time, TTmotion2.accel1);
% title('Head acceleration')
% ylabel('rad/s^2')
% 
% subplot(3,2,6)
% plot(TTmotion2.Time, TTmotion2.accel2);
% title('Cage acceleration')
% ylabel('rad/s^2')
% hold off;