%{
Converts (Aurora) data from euler angles to rotation data used
in a rotation matrix. Equations provided by Aurora
    
Nicole Devos for the WearME lab, Western University

ver 1.1
Nov 28, 2019

%}

function m = eulerConv(rvector)
% eulerConv   Converts (Aurora) data from euler angles to rotation data used
% in a rotation matrix.

% roll: Rz, Pitch: Ry, Yaw: Rx
% the rvector is Rz (1), Ry (2), Rx (3)
r = rvector(1); % was 3
p = rvector(2);
y = rvector(3); % was 1

m = zeros(3,3);
m(1,1) = cosd(r)*cosd(p);
m(1,2) = cosd(r)*sind(p)*sind(y) - sind(r)*cosd(y);
m(1,3) = cosd(r)*sind(p)*cosd(y) + sind(r)*sind(y);
m(2,1) = sind(r)*cosd(p);
m(2,2) = sind(r)*sind(p)*sind(y) + cosd(r)*cosd(y);
m(2,3) = sind(r)*sind(p)*cosd(y) - cosd(r)*sind(y);
m(3,1) = -sind(p);
m(3,2) = cosd(p)*sind(y);
m(3,3) = cosd(p)*cosd(y);

end
