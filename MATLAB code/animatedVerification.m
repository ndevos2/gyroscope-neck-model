%{

Animated position checker

Version 0.1
January 9th, 2022

%}

close all

% Exp to model
name = 'exp1-1';

% axes limit
XMAX = 250;
XMIN = -175;
YMIN = -200;
YMAX = 200;
ZMIN = -400;
ZMAX = 0; % -150;

% Data
load(strcat('MAT-Files/',name,'.mat'));
x = gnegate(TTa.s1Tx); % gnegate
y = TTa.s1Ty;
z = TTa.s1Tz;
x1 = gnegate(TTa.s2Tx); % gnegate
y1 = TTa.s2Ty;
z1 = TTa.s2Tz;

% initial sensor positions
das1 = mean(TTinit{1:50,4:6});
das2 = mean(TTinit{1:50,10:12});

% landmarks
load('MAT-Files/feb-13-distance-data.mat');
bx = -auroraToBase(1); % -
by = auroraToBase(2);
bz = auroraToBase(3);

sx = -auroraToSkullLandmark(1);
sy = auroraToSkullLandmark(2);
sz = auroraToSkullLandmark(3);

cx = -auroraToCageLandmark(1);
cy = auroraToCageLandmark(2);
cz = auroraToCageLandmark(3);

wx = -auroraToWheelLandmark(1);
wy = auroraToWheelLandmark(2);
wz = auroraToWheelLandmark(3);

% Create first frame of the animation
ani = scatter3(z(1),y(1),x(1),'*');
hold on;
axes = plot3([0;25],[0;0],[0;0],'-',[0;0],[0;-25],[0;0],'->',[0;0],[0;0],[0;100],'-^');
axes(1).MarkerFaceColor = axes(1).Color;
axes(2).MarkerFaceColor = axes(2).Color;
axes(3).MarkerFaceColor = axes(3).Color;
aur = scatter3(0,0,0,200,'filled','square');
ani2 = scatter3(z1(1),y1(1),x1(1),'or');
baxes = plot3([bz;bz+50],[by;by],[bx;bx],'->',[bz;bz],[by;by-25],[bx;bx],'->',[bz;bz],[by;by],[bx;bx+100],'-^');
baxes(1).MarkerFaceColor = baxes(1).Color;
baxes(2).MarkerFaceColor = baxes(2).Color;
baxes(3).MarkerFaceColor = baxes(3).Color;
base = scatter3(bz,by,bx,200,'filled','square');
skull = scatter3(sz,sy,sx,'x');
cage = scatter3(cz,cy,cx,'x');
wheel = scatter3(wz,wy,wx,'x');
pdas1 = scatter3(das1(3),das1(2),-das1(1),'b*');
pdas2 = scatter3(das2(3),das2(2),-das2(1),'or');

% Set the axes to fit the whole motion area
axis([ZMIN ZMAX YMIN YMAX XMIN XMAX])
xlabel('z')
ylabel('y')
zlabel('x')

% At each frame, add new position data point
for i = repmat(1:length(TTa.s1Tz), 1, 3)
    % Update the existing object's properties rather than creating a new one
    ani.XData = z(i);
    ani.YData = y(i);
    ani.ZData = x(i);
    ani2.XData = z1(i);
    ani2.YData = y1(i);
    ani2.ZData = x1(i);
    % Let you see the animation
    pause(0.01)
end

hold off;

%% Section for animating the tranformed part

figure;

% axes limits
XMAX = 400;
XMIN = -10;
YMIN = -10;
YMAX = 350;
ZMIN = -75;
ZMAX = 330; % -150;

% Transformation(s)
RA0 = Rz(pi)*Rx(pi/2); % Rotation matrix for Aurora (A) to base fram (0)
dA0 = auroraToBase; % position of Base Frame, measured with Aurora

% inverse (to convert Aurora data to Base frame)
H0A = [transpose(RA0), -transpose(RA0)*dA0; 0 0 0 1];

% TT of data to be transformed
TTtrans = removevars(TTa,{'s1Rx','s1Ry','s1Rz','s2Rx','s2Ry','s2Rz'});

% Location of Aurora
temp = H0A*[[1 0 0; 0 1 0; 0 0 1],[0;0;0];0,0,0,1];
aur0A = temp(1:3,4)';

% Location of skull landmark
temp = H0A*[[1 0 0; 0 1 0; 0 0 1],auroraToSkullLandmark;0,0,0,1];
skull0A = temp(1:3,4)';

% Location of cage landmark
temp = H0A*[[1 0 0; 0 1 0; 0 0 1],auroraToCageLandmark;0,0,0,1];
cage0A = temp(1:3,4)';

% Location of wheel landmark
temp = H0A*[[1 0 0; 0 1 0; 0 0 1],auroraToWheelLandmark;0,0,0,1];
wheel0A = temp(1:3,4)';

% starting positions of sensors
temp = H0A*[[1 0 0; 0 1 0; 0 0 1],mean(TTinit{1:50,4:6})';0,0,0,1];
idas1 = temp(1:3,4)';
temp = H0A*[[1 0 0; 0 1 0; 0 0 1],mean(TTinit{1:50,10:12})';0,0,0,1];
idas2 = temp(1:3,4)';

% Convert raw position data of both sensors to "wrt base"
for i=1:height(TTa)
    Ras1 = eulerConv(TTa{i,1:3});
    das1 = TTa{i,4:6}';
    Tas1 = [Ras1,das1;0 0 0 1];
    T02_1 = H0A*Tas1;
    
    Ras2 = eulerConv(TTa{i,7:9});
    das2 = TTa{i,10:12}';
    Tas2 = [Ras2 das2; 0 0 0 1];
    T02_2 = H0A*Tas2;

    TTtrans{i,:}=[T02_1(1:3,4)',T02_2(1:3,4)'];
end

% Set x, y, z coords for each sensor 
x = TTtrans.s1Tx;
y = TTtrans.s1Ty;
z = TTtrans.s1Tz;
x1 = TTtrans.s2Tx;
y1 = TTtrans.s2Ty;
z1 = TTtrans.s2Tz;

az = -aur0A(3);
ay = aur0A(2);
ax = aur0A(1);

% Create first frame of the animation
ani = scatter3(gnegate(z(1)),y(1),x(1),'*');
hold on;
axes = plot3([az;az+50],[ay;ay],[ax;ax],'->',[az;az],[ay;ay-25],[ax;ax],'->',[az;az],[ay;ay],[ax;ax+100],'-^');
axes(1).MarkerFaceColor = axes(1).Color;
axes(2).MarkerFaceColor = axes(2).Color;
axes(3).MarkerFaceColor = axes(3).Color;

aur = scatter3(az,ay,ax,200,'filled','square');

ani2 = scatter3(gnegate(z(1)),y(1),x(1),'or');

baxes = plot3([0;50],[0;0],[0;0],'->',[0;0],[0;-25],[0;0],'->',[0;0],[0;0],[0;100],'-^');
baxes(1).MarkerFaceColor = baxes(1).Color;
baxes(2).MarkerFaceColor = baxes(2).Color;
baxes(3).MarkerFaceColor = baxes(3).Color;

base = scatter3(0,0,0,200,'filled','square');
skull = scatter3(-skull0A(3),skull0A(2),skull0A(1),'x');
cage = scatter3(-cage0A(3),cage0A(2),cage0A(1),'x');
wheel = scatter3(-wheel0A(3),wheel0A(2),wheel0A(1),'x');

pidas1 = scatter3(-idas1(3),idas1(2),idas1(1),'b*');
pidas2 = scatter3(-idas2(3),idas2(2),idas2(1),'or');

% Set the axes to fit the whole motion area
axis([ZMIN ZMAX YMIN YMAX XMIN XMAX])
xlabel('-z')
ylabel('y')
zlabel('x')

% At each frame, add new position data point
for i = repmat(1:length(TTa.s1Tz), 1, 3)
    % Update the existing object's properties rather than creating a new one
    ani.XData = gnegate(z(i));
    ani.YData = y(i);
    ani.ZData = x(i);
    ani2.XData = gnegate(z1(i));
    ani2.YData = y1(i);
    ani2.ZData = x1(i);
    % Let you see the animation
    pause(0.01)
end

hold off;

%% angle animation

% Data
[anglex, angley] = pol2cart(deg2rad(convertedAngles(:,1)),200);
x = pos_s1(:,1);
y = pos_s1(:,2);

figure;

% Create first frame of the animation
ani = plot([0,x(1)],[0,y(1)]);
hold on;
ani2 = scatter(x(1), y(1), 100, '*');
ani3 = plot([0,anglex(1)],[0,angley(1)]);

%axis([200 -200 200 -200])

% At each frame, add new position data point
for i = repmat(1:length(convertedAngles), 1, 3)
    % Update the existing object's properties rather than creating a new one
    ani.XData = [0, x(i)];
    ani.YData = [0, y(i)];
    ani2.XData = x(i);
    ani2.YData = y(i);
    ani3.XData = [0, anglex(i)];
    ani3.YData = [0, angley(i)];
    % Let you see the animation
    pause(0.01)
end

hold off;

%% Check that the rate of change in the angles is constant

diff = zeros(length(convertedAngles),1);

for i=1:length(convertedAngles)
    sensorpos = rad2deg(cart2pol(x(i),y(i)));
    diff(i) = sensorpos - convertedAngles(i,1);
end

figure;
plot(diff);
minima = min(diff)
maxima = max(diff)