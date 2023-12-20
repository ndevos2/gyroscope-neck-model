function [sensor1,sensor2,chanAvgs] = initAvg(TTin)
%INITAVG Accepts a 12 column input timetable of Aurora data, averages every
%element in the position columns. Outputs the positions of sensors 1 and 2,
%optionally the initial average for each channel.
%   Uses the first 100 data points only. Detailed explanation goes here
% 
% Nicole Devos for the WearME Lab, Western University
% 
% Version 1.2
% December 14th, 2022

% Convert timetable into table, because MATLAB refuses to do 'mean' on a
% timetable for some dumb reason. Then you can use mean on the table.
Tin = timetable2table(TTin(1:100,:),'ConvertRowTimes',false);
avgs = mean(Tin{:,:});

% Set sensor averages based on Tx/y/z data
sensor1 = avgs(4:6).';
sensor2 = avgs(10:12).';

% Compile the channel averages in a READY TO PLOT vector.
index1 = 1:2:12;
index2 = 2:2:12;
chanAvgs(index1) = [avgs(4:6),flip(avgs(1:3))];
chanAvgs(index2) = [avgs(10:12),flip(avgs(7:9))];

end