%% DSP Project
% * Rachel Greenlee, Colten Humphrey
% * Digital Signal Processing
% * February 19, 2017
%% Clean up the environment
clear all;
format compact;
close all;
clc;

%% Start main
fs = 48000;
bits = 16;
channel = 2;
recObj = audiorecorder(fs,bits,channel);
set(recObj,'TimerPeriod',1,'TimerFcn',{@audioTimer});

recObj = audiorecorder;

for i=1:20
recordblocking(recObj, .5);
%# Store data in double-precision array.
sound = getaudiodata(recObj);    
    
timeValue = 0;
pitchValue = -0.4;
reverbValue = 0.3;
minfreq = 40;
maxfreq = 600;
ratios = 8^(timeValue-0.4*pitchValue); %the ratio of new frequency vs old frequency
n = 1024; %fft window size (default: 1024)
H = n/4; %fft window offset (default: n/4)
%fs = 22050;

%% Timescale
%Get the number of windows and the number of frequencies for the sfft
specgram = stft(sound,n,n,H);
[nFreq, nWindows] = size(specgram);

    %Put all the ratios in an array
    if(size(ratios,1) == 1)
        timeScale = 0:ratios:(nWindows-2);
    else
        cumul = 0;
        timeScale = 0;
        timeScaleInvRatio = 1;
        while(cumul < nWindows-2)
           timeScale = [timeScale cumul];
           ratio = ratios(ceil(cumul) + 1);
           cumul = ratio + cumul;
           invRatio = 1/ratio;
           timeScaleInvRatio = [timeScaleInvRatio invRatio];
        end
    end
    
    %Init the spectrogram with zeros
    outputSpectro = zeros(nFreq, length(timeScale));
    
    %Estimate phase offset for each window
    deltaPhaseAve = zeros(1, nFreq);
    deltaPhaseAve(2:nFreq) = (2 * pi * H) ./((nFreq*2)./(1:nFreq-1));
    
    %Init phase sum
    cumulPhase = angle(specgram(:,1));
    coloneCourante = 1;
    
    %Main loop for time scaling
    for t = timeScale 
        
       % All the columns
       interCols = specgram(:, floor(t) + [1 2]);
       
       %  position offset between the sample and the real value
       pos = t - floor(t);
       
       % interpolation of the 2 samples
       frame = (1-pos)*abs(interCols(:,1)) + pos*abs(interCols(:,2));
       
       % Phase difference between observed value and estimated value
       deltaPhase = angle(interCols(:,2)) - angle(interCols(:,1)) - deltaPhaseAve';
       
       %Scale the phase to a number between -pi and pi
       deltaPhase = deltaPhase - 2 * pi * round(deltaPhase/(2*pi));
       
       %Fix the phase
       outputSpectro(:, coloneCourante) = frame .* exp(1i*cumulPhase);
       
       %Next iteration
       cumulPhase = cumulPhase + deltaPhaseAve' + deltaPhase;
       coloneCourante = coloneCourante + 1;
    end
    
    %output signal
    sound = istft(outputSpectro,n,n,H)';

%% Pitch shift
[N,D] = rat(8^(-0.4*pitchValue));
sound = resample(sound,N,D);

%% Reverb
	delaytime = 120;
	M=floor(delaytime./1000.*fs);

	b = zeros(M,1);
	b(1) = 1;
	b(M) = -reverbValue;
	a = 1;

	h = impz(b,a);
	sound = conv(sound,h);

%% Stft
    sound = stft(sound, 1024, 1024, 256, fs);
%%
    display(sound);
    sound(sound, fs, bits);

end

%# Plot the waveform.
plot(player);

