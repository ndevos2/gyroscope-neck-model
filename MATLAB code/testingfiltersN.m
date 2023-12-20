function testingfiltersN(data)

%Trying to come up with the right filter to use for the position and force
%data.
% Modified from testingfilters.m by Ana Luisa Trejos, Feb, 2012. Modifed by
% Nicole Devos, May 18, 2021.

%assign variables and zero from starting position
dataX = data(:,4)-data(1,4);
dataZ = data(:,6)-data(1,6);
s = 0.025; %this is my sampling rate for position
Fs = 40; %this is the sampling frequency for position
%s = 0.002; %this is my sampling rate for force
%Fs = 500; %this is the sampling frequency for force

l = length(data);
t = l*s; % convert count to seconds - this is the duration

%plot the position
tt = linspace(0,t,l);
figure
subplot(2,1,1), plot(tt, dataX)
title('Tx')
subplot(2,1,2), plot(tt, dataZ)
title('Tz')
sgtitle('Original position signals')

%figure out cutoff frequency
NFFT = 2^nextpow2(l); % Next power of 2 from length of y
Y = fft(dataX,NFFT)/l;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')


% Power spectrum plot
fftx = fft(dataX,NFFT);
NumUniquePts = ceil((NFFT+1)/2);
fftx = fftx(1:NumUniquePts);
mx = abs(fftx);
mx = mx/length(dataX);
mx = mx.^2; 
if rem(NFFT, 2) % odd nfft excludes Nyquist point 
  mx(2:end) = mx(2:end)*2;
else
  mx(2:end -1) = mx(2:end -1)*2;
end
f = (0:NumUniquePts-1)*Fs/NFFT;
% Generate the plot, title and labels.
figure
plot(f,10*log10(mx));
title('Power Spectrum');
xlabel('Frequency (Hz)'); 
ylabel('Power (dB)');


%filter
%[b,a] = butter(4,0.005); %for force
% [b,a] = butter(4,0.025); %for position - 0.5 Hz
[b,a] = butter(4,0.05); % 1 Hz - lowpass/(fs/2)
X = filtfilt(b,a,dataX);
Y = filtfilt(b,a,dataZ);
% X = fastsmooth(dataX,w,2,1);
% Y = fastsmooth(dataY,w,2,1);
% Z = fastsmooth(dataZ,w,2,1);

%plot the position
tt = linspace(0,t,l);
% figure
% subplot(3,1,1), plot(tt, dataX, 'r', tt, X, 'b')
% title('Filtered signal vs. original signal')
% subplot(3,1,2), plot(tt, dataY, 'r', tt, Y, 'b')
% figure
% plot(tt, dataZ, 'r', tt, Z, 'b')
% axis([110 140 20 150])


% Power spectrum plot
fftx = fft(X,NFFT);
NumUniquePts = ceil((NFFT+1)/2);
fftx = fftx(1:NumUniquePts);
mx = abs(fftx);
mx = mx/length(X);
mx = mx.^2; 
if rem(NFFT, 2) % odd nfft excludes Nyquist point 
  mx(2:end) = mx(2:end)*2;
else
  mx(2:end -1) = mx(2:end -1)*2;
end
f = (0:NumUniquePts-1)*Fs/NFFT;
% Generate the plot, title and labels.
hold on
plot(f,10*log10(mx),'r');
title('4th Order Butterworth - Power Spectrum Plot')

%first derivative
DX=diff(X)/s;
DY=diff(Y)/s;

figure
subplot(2,1,1), plot(tt(1:end-1), DX)
title('DX')
sgtitle('4th Order Butterworth - First derivative')
subplot(2,1,2), plot(tt(1:end-1), DY)
title('DZ')

% figure
% plot(tt(1:end-1),VelA1)

%second derivative
DdX=diff(DX)/s;
DdY=diff(DY)/s;

figure
subplot(2,1,1), plot(tt(1:end-2), DdX)
title('DdX')
sgtitle('4th Order Butterworth - Second derivative')
subplot(2,1,2), plot(tt(1:end-2), DdY)
title('DdZ')

% Power spectrum plot
l = length(DX);
NFFT = 2^nextpow2(l);
fftx = fft(DX,NFFT);
NumUniquePts = ceil((NFFT+1)/2);
fftx = fftx(1:NumUniquePts);
mx = abs(fftx);
mx = mx/length(DX);
mx = mx.^2; 
if rem(NFFT, 2) % odd nfft excludes Nyquist point 
  mx(2:end) = mx(2:end)*2;
else
  mx(2:end -1) = mx(2:end -1)*2;
end
f = (0:NumUniquePts-1)*Fs/NFFT;
% Generate the plot, title and labels.
hold on
figure
plot(f,10*log10(mx),'r');
title('4th Order Butterworth - Power Spectrum Plot - First Derivative')


%% repeat for 6th order filter

% Power spectrum plot
fftx = fft(dataX,NFFT);
NumUniquePts = ceil((NFFT+1)/2);
fftx = fftx(1:NumUniquePts);
mx = abs(fftx);
mx = mx/length(dataX);
mx = mx.^2; 
if rem(NFFT, 2) % odd nfft excludes Nyquist point 
  mx(2:end) = mx(2:end)*2;
else
  mx(2:end -1) = mx(2:end -1)*2;
end
f = (0:NumUniquePts-1)*Fs/NFFT;
% Generate the plot, title and labels.
figure
plot(f,10*log10(mx));
title('Power Spectrum');
xlabel('Frequency (Hz)'); 
ylabel('Power (dB)');

% [b,a] = butter(6,0.025); %for position - [Hz]/(40/2) - 0.5 Hz
[b,a] = butter(6,0.05); %for position - [Hz]/(40/2) - 1 Hz
X = filtfilt(b,a,dataX);
Y = filtfilt(b,a,dataZ);

l = length(data);
tt = linspace(0,t,l);

% Power spectrum plot
fftx = fft(X,NFFT);
NumUniquePts = ceil((NFFT+1)/2);
fftx = fftx(1:NumUniquePts);
mx = abs(fftx);
mx = mx/length(X);
mx = mx.^2; 
if rem(NFFT, 2) % odd nfft excludes Nyquist point 
  mx(2:end) = mx(2:end)*2;
else
  mx(2:end -1) = mx(2:end -1)*2;
end
f = (0:NumUniquePts-1)*Fs/NFFT;
% Generate the plot, title and labels.
hold on
plot(f,10*log10(mx),'r');
title('6th Order Butterworth - Power Spectrum Plot')

%first derivative
DX=diff(X)/s;
DY=diff(Y)/s;

figure
subplot(2,1,1), plot(tt(1:end-1), DX)
title('DX')
sgtitle('6th Order Butterworth - First derivative')
subplot(2,1,2), plot(tt(1:end-1), DY)
title('DZ')

% figure
% plot(tt(1:end-1),VelA1)

%second derivative
DdX=diff(DX)/s;
DdY=diff(DY)/s;

figure
subplot(2,1,1), plot(tt(1:end-2), DdX)
title('DdX')
sgtitle('6th Order Butterworth - Second derivative')
subplot(2,1,2), plot(tt(1:end-2), DdY)
title('DdZ')

% Power spectrum plot
l = length(DX);
NFFT = 2^nextpow2(l);
fftx = fft(DX,NFFT);
NumUniquePts = ceil((NFFT+1)/2);
fftx = fftx(1:NumUniquePts);
mx = abs(fftx);
mx = mx/length(DX);
mx = mx.^2; 
if rem(NFFT, 2) % odd nfft excludes Nyquist point 
  mx(2:end) = mx(2:end)*2;
else
  mx(2:end -1) = mx(2:end -1)*2;
end
f = (0:NumUniquePts-1)*Fs/NFFT;
% Generate the plot, title and labels.
hold on
figure
plot(f,10*log10(mx),'r');
title('6th Order Butterworth - Power Spectrum Plot - First Derivative')

