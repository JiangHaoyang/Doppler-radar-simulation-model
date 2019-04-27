clear all;
close all;
% switches of different modules
playAnimation='y';
runSpectrumMethod='n';
runDopplerDelayMethod='n';generateRetuenSignal='y';generateFrames='y';
% if runDopplerDelayMethod==n, then 'generateRetuenSignal' and 'generateFrames'
% do not work.
%% set parameters
% car parameters
carLength=4.5; % length in m
carWidth=2; % width in m
radius=0.4; % wheel radius in m
thickness=0.1; % thickness of ellipsoids
ringsInWheel=2; 
PointsInRing=5;
carParameters=[carLength,carWidth,radius,thickness,ringsInWheel,PointsInRing];

% trajectory
trajectory=@(t) 3*(t/2.5-sin(t/2.5));% the car runs on the x axis.

% signal parameters
sample_rate=1e10; % samplin Frequency in Hz
startF = 4e9; % start frequency in Hz
endF=4.5e9; % stop frequency in Hz
carrierFrequency = (startF+endF)/2;
BandWidth=endF-startF; % bandwidth in Hz
c=3e8; % speed of light in m/s
Signal_duration=1.5e-7; % pulse repetition period in s
Pulse_duration=0.2*Signal_duration; % pulse width in s
Pulse_time=0:1/sample_rate:Signal_duration; % short time covering one PRP

% Signal=@(t) 1/sqrt(Pulse_duration)...
%     *RectangleWave(t,Pulse_duration)...
%     .*exp(1i*pi*BandWidth/Pulse_duration*t.^2)...
%     .*exp(1i*startF*2*pi*t);
Signal=@(t) RectangleWave(t,Pulse_duration)...
    .*exp(1i*pi*BandWidth/Pulse_duration*t.^2)...
    .*exp(1i*startF*2*pi*t);

% time parameters - long time
startTime=0;
endTime=16;
LongSampleRate=500;
longTimeParameters=[startTime,endTime,LongSampleRate];
dT=1/LongSampleRate;


% radar location
radarLocation=[9,7,7];

% calculate position data
wheelPositionData = carModel( carParameters,longTimeParameters,trajectory );

%% play the animation
if playAnimation=='y'
    AXIS=[-7,23,-3,8,-1,9.5];
    [ animation ] = Car_animation( wheelPositionData,carParameters,radarLocation,AXIS );
end

%% spetrum method
if runSpectrumMethod=='y'
    [ time,doppler,spectrogram1 ] = ChenMethod( carrierFrequency,longTimeParameters,radarLocation,carParameters,wheelPositionData );
end

%% doppler-delay method
if runDopplerDelayMethod=='y'
    shortTimeParameters=[sample_rate,startF,endF,Signal_duration,Pulse_duration,c];
    controlFlag=[generateRetuenSignal,generateFrames];
    [DopplerDelayAnimation,spectrogram2,spectrogram_timeAxis,spectrogram_dopplerAxis,spectrogramCmplx]=doppler_delayMethod...
        ( Signal,shortTimeParameters,longTimeParameters,...
        carParameters,wheelPositionData,radarLocation,controlFlag);
    if (spectrogram2~=0)
        fig=figure('name','spectrum2');
        imagesc(spectrogram_timeAxis,spectrogram_dopplerAxis,fftshift(20*log10(spectrogram2)-max(max(20*log10(spectrogram2))),1));
        colormap('jet');
        colorbar;
        caxis([-45,0]);
        xlabel('time(s)');
        ylabel('doppler(Hz)');
        axis xy;
    end
end