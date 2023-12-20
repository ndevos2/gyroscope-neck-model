function R = Rz(rotAngle)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rotAngle = sym(rotAngle);

R = [cos(rotAngle) -sin(rotAngle) 0; sin(rotAngle) cos(rotAngle) 0; 0 0 1];

R = double(R);

end