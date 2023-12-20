function [PVinF,PLinF,TVinF,TLinF] = multiPT(TTin,Name)
%MULTIPT graphs each channel in the dataset with the local minima and
%maxima indicated and numbered. Returns the location and values for each
%peak/trough in each channel.
%   The function assumes the TT input is 2-sensor Aurora data if the TT
%   has 12 variables, or the TT input is ATI force/torque data if the TT
%   has 6 variables. Using the findpeaks formula, calculates the Peak
%   Values and Peak Locations of each channel of the signal. Using the
%   gnegate formula (general negation) and findpeaks, calculated the Trough
%   Values and Trough Locations of the signal. Plots each channel findpeaks
%   as a subplot, and plots the trough values on top.
%   The second input, Name, is used to name the figure window and title
%   the graph.

if width(TTin)==12
    channel = {'s1Tx','s2Tx','s1Ty','s2Ty','s1Tz','s2Tz','s1Rx','s2Rx','s1Ry','s2Ry','s1Rz','s2Rz'};

    figure('Name',Name);
    sgtitle(Name);

    for i = 1:length(channel)
        subplot(6,2,i);
        [PVinF.(channel{i}),PLinF.(channel{i})] = findpeaks(TTin.(channel{i}));
        [TVinF.(channel{i}),TLinF.(channel{i})] = findpeaks(gnegate(TTin.(channel{i})));
        findpeaks(TTin.(channel{i}))
        text(PLinF.(channel{i})+12,PVinF.(channel{i})+.5,num2str((1:numel(PVinF.(channel{i})))'))
        text(TLinF.(channel{i})+12,gnegate(TVinF.(channel{i}))-.5,num2str((1:numel(TVinF.(channel{i})))'))
        hold on;
        scatter(TLinF.(channel{i}),gnegate(TVinF.(channel{i})),'^',"filled")
        title(channel{i})
        hold on;
    end
end

if width(TTin)==6
    channel = {'Fx','Mx','Fy','My','Fz','Mz'};

    figure('Name',Name);
    sgtitle(Name);

    for i = 1:length(channel)
        subplot(3,2,i);
        [PVinF.(channel{i}),PLinF.(channel{i})] = findpeaks(TTin.(channel{i}));
        [TVinF.(channel{i}),TLinF.(channel{i})] = findpeaks(gnegate(TTin.(channel{i})));
        findpeaks(TTin.(channel{i}))
        text(PLinF.(channel{i})+0,PVinF.(channel{i})+0,num2str((1:numel(PVinF.(channel{i})))'))
        text(TLinF.(channel{i})+0,gnegate(TVinF.(channel{i}))-0,num2str((1:numel(TVinF.(channel{i})))'))
        hold on;
        scatter(TLinF.(channel{i}),gnegate(TVinF.(channel{i})),'^',"filled")
        title(channel{i})
        hold on;
    end
end

end