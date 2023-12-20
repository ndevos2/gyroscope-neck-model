function error = checkError(TTin)
%CHECKERROR Checks for errors in the TT input; presented as a percent of
%total data points
%   Checks for percent errors in the TT input (indicated by a number that
%   is smaller than -1e10, which is the same threshold for the fixerrors
%   function). Presented as a percentage of total data points. Calculates
%   for each channel, either 6 or 12 (depending on the number of sensors in
%   the dataset).

total = height(TTin);

error.s1Tx = sum(TTin.s1Tx<-1*10^10)/total;
error.s1Ty = sum(TTin.s1Ty<-1*10^10)/total;
error.s1Tz = sum(TTin.s1Tz<-1*10^10)/total;
error.s1Rx = sum(TTin.s1Rx<-1*10^10)/total;
error.s1Ry = sum(TTin.s1Ry<-1*10^10)/total;
error.s1Rz = sum(TTin.s1Rz<-1*10^10)/total;

if width(TTin) == 12
   error.s2Tx = sum(TTin.s2Tx<-1*10^10)/total;
   error.s2Ty = sum(TTin.s2Ty<-1*10^10)/total;
   error.s2Tz = sum(TTin.s2Tz<-1*10^10)/total;
   error.s2Rx = sum(TTin.s2Rx<-1*10^10)/total;
   error.s2Ry = sum(TTin.s2Ry<-1*10^10)/total;
   error.s2Rz = sum(TTin.s2Rz<-1*10^10)/total;
end
end