clear all;
% close all;
clc;
% switches of different modules
playAnimation='y';
runChenMethod='n';
runDopplerDelayMethod='n';generateRetuenSignal='y';generateFrames='y';
% if runDopplerDelayMethod==n, then 'generateRetuenSignal' and 'generateFrames'
% do not work.
%% wind turbine parameters
tower1_height=100;tower1_radius=4.5;
tower2_height=30 ;tower2_radius=4;
tower3_height=20 ;tower3_radius=3.5;
nacelle_x=20; nacelle_y=7.5; nacelle_z=8;
nacelle_bias_x=5;nacelle_bias_y=nacelle_y/2;
rotor_radius=4;rotorCneterBeyondTheNacelle=3;
blade_1_height=20;blade_1_radius=2.5;
blade_2_x=20;% 20(typical value)
blade_2_y=5;% 5
blade_2_z=1;% 1
blade_2_rotate=20/180*pi;%20 degree
blade_3_x=20;% 20
blade_3_y=4;% 4
blade_3_z=0.5;%0.5
blade_3_rotate=40/180*pi;%40 degree
blade_4_x=20;% 20
blade_4_y=2;% 2
blade_4_z=0.3;% 0.3
blade_4_rotate=60/180*pi;%60 degree
RotateSpeed=10;% rpm (rounds per min)

ModelParameters=[tower1_height , tower1_radius  , 0 , 0 ;        % bottom of the tower
                 tower2_height , tower2_radius  , 0 , 0 ;        % middle of the tower
                 tower3_height , tower3_radius  , 0 , 0 ;        % top of the tower
                 nacelle_x     , nacelle_y,     nacelle_z ,0 ;  % size of the nacelle
                 nacelle_bias_x, nacelle_bias_y , 0 , 0 ;        % tower Top At The Bottom Of TheNacelle
                 rotor_radius  , rotorCneterBeyondTheNacelle,0 ,0 ;
                 blade_1_height, blade_1_radius ,0 , 0 ;
                 blade_2_x     , blade_2_y      ,blade_2_z , blade_2_rotate;
                 blade_3_x     , blade_3_y      ,blade_3_z , blade_3_rotate;
                 blade_4_x     , blade_4_y      ,blade_4_z , blade_4_rotate;
                 RotateSpeed   , 0 , 0 , 0 ];
ModelParameters(1:7,:)=ModelParameters(1:7,:)/3;
ModelParameters(8:10,1:3)=ModelParameters(8:10,1:3)/3;
%% signal parameters
sample_rate=2.5e9; % samplin Frequency in Hz
startF = 4e9; % start frequency in Hz
endF=4.5e9; % stop frequency in Hz
carrierFrequency = (startF+endF)/2;
BandWidth=endF-startF; % bandwidth in Hz
c=3e8; % speed of light in m/s
Signal_duration=7e-7; % pulse repetition period in s
Pulse_duration=0.2*Signal_duration; % pulse width in s
Pulse_time=0:1/sample_rate:Signal_duration; % short time covering one PRP

% Signal=@(t) 1/sqrt(Pulse_duration)...
%     *RectangleWave(t,Pulse_duration)...
%     .*exp(1i*pi*BandWidth/Pulse_duration*t.^2)...
%     .*exp(1i*startF*2*pi*t);
Signal=@(t) RectangleWave(t,Pulse_duration)...
    .*exp(1i*pi*BandWidth/Pulse_duration*t.^2)...
    .*exp(1i*startF*2*pi*t);

%% long time parameters
startTime=0;
endTime=7;
LongSampleRate=1600;
longTimeParameters=[startTime,endTime,LongSampleRate];
dT=1/LongSampleRate;

%% radar location
radarLocation=[10,-20,2];

%% calculate position data
windTurbine_modelData = WindTurbineModel( ModelParameters,longTimeParameters );
%% play the animation
if playAnimation=='y'
    AXIS=[-12,15,-30,30,-5,85];
    [ animation ] = WindTurbine_animation( windTurbine_modelData,radarLocation,AXIS );
end
%% Chen method
if runChenMethod=='y'
    [ spectrogram,doppler,time ] = ChenMethod_forWindTurbine( windTurbine_modelData,radarLocation,longTimeParameters,carrierFrequency );
end

%% doppler-delay method
if runDopplerDelayMethod=='y'
    shortTimeParameters=[sample_rate,startF,endF,Signal_duration,Pulse_duration,c];
    controlFlag=[generateRetuenSignal,generateFrames];
    [DopplerDelayAnimation,spectrogram2,spectrogram_timeAxis,spectrogram_dopplerAxis,spectrogramCmplx]=Doppler_delayMethod_forWindTurbine...
        ( Signal,shortTimeParameters,longTimeParameters,...
        windTurbine_modelData,radarLocation,controlFlag);
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