function [TToutput] = fixerrors(TTinput)
%fixerrors Find errors in timetable (e28) and fixes them
%   Locates all instances of - e10, replaces with NaN, and uses "fill
%   missing" to interpolate the missing data points. Must be in matrix
%   format vs timetable format; takes out the variable info and puts back
%   into TT format before returning TToutput.
%
% Nicole Devos for the WearME lab, Western University
% Ver 1.1, July 6, 2022
%

%input = [TTinput.s1Rz,TTinput.s1Ry,TTinput.s1Rx,TTinput.s1Tx,TTinput.s1Ty,TTinput.s1Tz,TTinput.s2Rz,TTinput.s2Ry,TTinput.s2Rx,TTinput.s2Tx,TTinput.s2Ty,TTinput.s2Tz];
input = TTinput.Variables;
input(input<-1*10^10) = NaN;
input = fillmissing(input,'linear');

% Aurora frequency is 40
%TToutput = array2timetable(input,'SampleRate',40);
%TToutput.Properties.VariableNames = ["s1Rz","s1Ry","s1Rx","s1Tx","s1Ty","s1Tz","s2Rz","s2Ry","s2Rx","s2Tx","s2Ty","s2Tz"];
TTinput{:,:} = input;
TToutput = TTinput;

end

