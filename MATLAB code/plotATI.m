function [f1] = plotATI(TTin,Name,TTsecond)
%PLOTATI Plots the individual force and torque data from a Timetable input
%   Plots the individual force and torque data of an ATI Timetable.
%   Force and Torque side-by-side for easier comparison.
%   Can allow for signals to be overlaid for comparison. For example, a
%   filtered and raw signal.
%   
%   Nicole Devos for the WearME Lab, Western University
%   Version 1.0, December 14, 2022

arguments
    TTin {mustBeA(TTin,"timetable")}
    Name 
    TTsecond = 0
end

channels1 = {"Fx","Mx","Fy","My","Fz","Mz"};

f1 = figure('Name',Name);
sgtitle(Name);
for i = 1:6
    subplot(3,2,i);
    plot(TTin.Time,TTin.(channels1{i}))
    title(channels1{i})
    hold on;
end
if istimetable(TTsecond)
    for i = 1:6
        subplot(3,2,i);
        plot(TTsecond.Time,TTsecond.(channels1{i}))
    end
end