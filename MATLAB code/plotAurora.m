function [f1] = plotAurora(TTin,Name,TTsecond)
%PLOTAURORA Plots the individual position and rotation data from a
%Timetable input
%   Plots the individual position and rotation data of an Aurora Timetable.
%   Sensor 1 and 2 side-by-side for comparison.
%   Can allow for two signals to be overlaid - for example, to compare a
%   raw and filtered signal.
%   
%   Nicole Devos for the WearME Lab, Western University
%   Version 1.3, December 09, 2022

arguments
    TTin {mustBeA(TTin,"timetable")}
    Name 
    TTsecond = 0
end

numSensors = width(TTin);
channels1 = {"s1Tx","s1Ty","s1Tz","s1Rx","s1Ry","s1Rz"};
channels2 = {"s1Tx","s2Tx","s1Ty","s2Ty","s1Tz","s2Tz","s1Rx","s2Rx","s1Ry","s2Ry","s1Rz","s2Rz"};

if numSensors == 6
    f1 = figure('Name',Name);
    sgtitle(Name);
    for i = 1:numSensors
        subplot(6,1,i);
        plot(TTin.Time,TTin.(channels1{i}))
        title(channels1{i})
        hold on;
    end
    if istimetable(TTsecond)
        for i = 1:numsensors
            subplot(6,1,i);
            plot(TTsecond.Time,TTsecond.(channels1{i}))
        end
    end


else
    f1 = figure('Name',Name);
    sgtitle(Name);
    for i = 1:numSensors
        subplot(6,2,i);
        plot(TTin.Time,TTin.(channels2{i}))
        title(channels2{i})
        hold on;
    end
    if istimetable(TTsecond)
        for i = 1:numSensors
            subplot(6,2,i)
            plot(TTsecond.Time,TTsecond.(channels2{i}))
        end
    elseif TTsecond ~= 0
        for i = 1:numSensors
            subplot(6,2,i)
            yline(TTsecond(i),'-r')
        end
    end

end

end

