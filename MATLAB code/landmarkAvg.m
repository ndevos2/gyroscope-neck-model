function [position] = landmarkAvg(filename)
%LANDMARKAVG This function averages every entry in the D:F columns and
%reformats to be a position vector
%   This function takes a filename from a CLEANED aurora file and
%imports the matrix and averages each column, and reformats into a position
%vector [i; j; k;]

% Confirm that MATLAB is reading the time column as the A column
file = readmatrix(filename,'Range','E:G');
position = mean(file(2:end,:)).';

end

