function [ spectrogram,doppler,time ] = ChenMethod_forWindTurbine( model_data,radarLocation,timeParameters,carrierFrequency )
%timeParameters=[startTime,endTime,sampleRate];
% model_data is a struct contains information of tower, nacelle, rotor and
% blades.
%
% *****
% model_data.tower is a 3*3*3 matrix
% model_data.tower(:,:,k) contains information of the k(th) part fo of the
% tower.
% Example: model_data.tower(:,:,k)=[P1;P2;[height,radius,0]]; P1 and P2 are
% in row vector.
%
% *****
% model_data.nacelle=[nacelleCenter;nacelleSize]
% size(model_data.nacelle) is [2,3];
%
% *****
% model_data.rotor=[rotorCenter;[radius,0,0]]
% size(model_data.rotor) is [2,3]
%
% *****
% model_data.bladeInfo=ModelParameters(7:10,:)
%
% *****
% one blade is modeled by 4 parts:
% P1---part1---P2---part2---P3---part3---P4---part4---P5
%
% *****
% model_data.bladePosition contain all information of P1 to P5 during the
% time.
%
% model_data.bladePosition is a <3*length(time)*3*5> 4 dimensions matrix
% the third number indicates the blade number(we have 3 blades in this model)
% the forth number indicates the P(we have 5 blades in this model)
% the first number indicates the axis(x,y,z)
% the second number indecates the time
% Example:
% model_data.bladePosition(:,(i),2,3)=[7,8,9]
% which means that the coordinate of the P3 of the 2nd blade at time(i) is [7,8,9]
%% function execuation

% set radar parameters
c=3e8;% light speed in m/s
Signal=@(t) exp(1i*carrierFrequency*2*pi*t);
radar_x=radarLocation(1);
radar_y=radarLocation(2);
radar_z=radarLocation(3);

% set time
startTime=timeParameters(1);
endTime=timeParameters(2);
sampleRate=timeParameters(3);
dT=1/sampleRate;
time=startTime:dT:endTime;

% generate raw data
returned_signal=zeros(size(time));

%% calculate signal returned from the tower
for i=1:3
    P1=model_data.tower(1,:,i);
    P2=model_data.tower(2,:,i);
    P=(P1+P2)/2;
    dx=P(1)-radar_x;
    dy=P(2)-radar_y;
    dz=P(3)-radar_z;
    distance=sqrt(dx.^2+dy.^2+dz.^2);
    rcs=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% not write
    amp=sqrt(rcs);
    returned_signal=returned_signal+amp.*Signal(time-2*distance/c);
end
%% calculate signal returned from the nacelle
P=model_data.nacelle(1,:);
dx=P(1)-radar_x;
dy=P(2)-radar_y;
dz=P(3)-radar_z;
distance=sqrt(dx.^2+dy.^2+dz.^2);
rcs=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% not write
amp=sqrt(rcs);
returned_signal=returned_signal+amp.*Signal(time-2*distance/c);
%% calculate signal returned from the roter
P=model_data.rotor(1,:);
dx=P(1)-radar_x;
dy=P(2)-radar_y;
dz=P(3)-radar_z;
distance=sqrt(dx.^2+dy.^2+dz.^2);
rcs=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% not write
amp=sqrt(rcs);
returned_signal=returned_signal+amp.*Signal(time-2*distance/c);
%% calculate signal returned from the blades
for m=1:3
    for n=1:4
        P1=model_data.bladePosition(:,:,m,n);
        P2=model_data.bladePosition(:,:,m,n+1);
        P=(P1+P2)/2;
        dx=P(1,:)-radar_x;
        dy=P(2,:)-radar_y;
        dz=P(3,:)-radar_z;
        distance=sqrt(dx.^2+dy.^2+dz.^2);
        rcs=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% not write
        amp=sqrt(rcs);
        returned_signal=returned_signal+amp.*Signal(time-2*distance/c);
    end
end

distance_ref=sqrt(radar_x.^2+radar_y.^2+radar_z.^2);
ref_signal=Signal(time-2*distance_ref/c);
% cancel the information of the carry frequency
signal=returned_signal./ref_signal;

%% code copy from the Chen.V's book

% micro-Doppler signature
f = signal;
np = length(time);
F = 1/dT;
T=endTime;
dF=1/T;
wd = 512;
wdd2 = wd/2;
wdd8 = wd/8;
ns = np/wd;

% calculate time-frequency micro-Doppler signature
disp('Calculating segments of TF distribution ...')
for k = 1:ns
    disp(strcat('  segment progress: ',num2str(k),'/',num2str(round(ns))))
    sig(1:wd,1) = f(1,(k-1)*wd+1:(k-1)*wd+wd);
    TMP = stft(sig,16);
    TF2(:,(k-1)*wdd8+1:(k-1)*wdd8+wdd8) = TMP(:,1:8:wd);
end
TF = TF2;
disp('Calculating shifted segments of TF distribution ...')
TF1 = zeros(size(TF));
for k = 1:ns-1
    disp(strcat('  shift progress: ',num2str(k),'/',num2str(round(ns-1))))
    sig(1:wd,1) = f(1,(k-1)*wd+1+wdd2:(k-1)*wd+wd+wdd2);
    TMP = stft(sig,16);
    TF1(:,(k-1)*wdd8+1:(k-1)*wdd8+wdd8) = TMP(:,1:8:wd);
end
disp('Removing edge effects ...')
for k = 1:ns-1
    TF(:,k*wdd8-8:k*wdd8+8) = ...
        TF1(:,(k-1)*wdd8+wdd8/2-8:(k-1)*wdd8+wdd8/2+8);
end

% display final time-frequency signature
fig=figure('name','spectrum1');
colormap(jet)
imagesc([0,T],[-F/2,F/2],20*log10(fftshift(abs(TF),1)+eps))
xlabel('Time (s)')
ylabel('Doppler (Hz)')
title('Micro-Doppler Signature of Rotating Rotor Blades')
axis xy
clim = get(gca,'CLim');
set(gca,'CLim',clim(2) + [-50 0]);
colorbar
drawnow

%% return parameters
%time=time;
doppler=linspace(-F/2,F/2,length(TF(:,1)));
spectrogram=TF;
end

