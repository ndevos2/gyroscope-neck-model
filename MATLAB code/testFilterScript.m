% load data

auroraFreq = 40;
% name = 'FEB 13 INTERP DATA/aurora-interp-feb-13-zero-degree-motor-test.csv';
% name = 'FEB 13 CLIPPED DATA/aurora-front-feb-13-exp1-1.csv';
name = 'FEB 13 INTERP DATA/aurora-interp-feb-13-side-side-test.csv';


TT = readtimetable(name);

testingfiltersN(TT{1:end,1:end})