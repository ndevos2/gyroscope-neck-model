function R = Ry(rotAngle)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

rotAngle = sym(rotAngle);

R = [cos(rotAngle) 0 sin(rotAngle); 0 1 0; -sin(rotAngle) 0 cos(rotAngle)];

R = double(R);

end