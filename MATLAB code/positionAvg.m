function [position] = positionAvg(filename,column)
%POSITIONAVG This function averages every entry in the specified columns and
%reformats to be a position vector
%   This function takes a filename from a CLEANED aurora file and
%imports the matrix and averages each column, and reformats into a position
%vector [i; j; k;]

file = readmatrix(filename,'Range',column);
position = mean(file(2:end,:));

end
