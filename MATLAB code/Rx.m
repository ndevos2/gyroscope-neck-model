function R = Rx(rotAngle)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rotAngle = sym(rotAngle);

R = [1 0 0; 0 cos(rotAngle) -sin(rotAngle); 0 sin(rotAngle) cos(rotAngle)];

R = double(R);

end