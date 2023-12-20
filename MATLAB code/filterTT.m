function TTout = filterTT(TTin,freq,lowhz)
%FILTERTT uses a butterworth filter on the input timetable
%   Uses the provided low pass frequency for a butterworth filter (lowhz),
%   and the input frequency of the TT (freq). Note that  the Aurora
%   frequency = 40Hz, ATI= 1000/16 Hz
% 
% Nicole Devos for the WearME lab, Western University
% Ver 1.1, December 8, 2022

fs = freq;
LP = lowhz/(fs/2);
[b,a] = butter(6,LP,'low');

TTout = TTin;

TTout{:,:} = filtfilt(b,a,TTin{:,:});

end