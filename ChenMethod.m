function [ time,doppler,spectrogram ] = ChenMethod( carrierFrequency,timeParameters,radarLocation,carParameters,wheelPositionData )
%timeParameters=[startTime,endTime,sampleRate];
%carParameters=[carLength,carWidth,radius,thickness,ringsInWheel,pointsInRing];

% set radar parameters
c=3e8;% light speed in m/s
Signal=@(t) exp(1i*carrierFrequency*2*pi*t);
radar_x=radarLocation(1);
radar_y=radarLocation(2);
radar_z=radarLocation(3);

% set the wheel paramaters;
% radius=carParameters(3);
thickness=carParameters(4);
rings=carParameters(5);
points=carParameters(6);

% set time
startTime=timeParameters(1);
endTime=timeParameters(2);
sampleRate=timeParameters(3);
dT=1/sampleRate;
time=startTime:dT:endTime;

% generate raw data
returned_signal=zeros(size(time));
%% calculate signal returned from axles
for j=1:2
    for k=1:2
        index=(j-1)*2+k;
        P1=wheelPositionData(index).center;
        if index==4
            P2=wheelPositionData(1).center;
        else
            P2=wheelPositionData(index+1).center;
        end
        P=(P1+P2)/2;
        dx=P(1,:)-radar_x;
        dy=P(2,:)-radar_y;
        dz=P(3,:)-radar_z;
        distance=sqrt(dx.^2+dy.^2+dz.^2);
        % calculate the RCS
        rcs=zeros(size(time));
        for i=1:length(time)
            A=[dx(i),dy(i),dz(i)];
            B=[P1(1,i)-P2(1,i),P1(2,i)-P2(2,i),P1(3,i)-P2(3,i)];
            A_dot_B = dot(A,B);
            A_sum_sqrt = sqrt(sum(A.*A));
            B_sum_sqrt = sqrt(sum(B.*B));
            ThetaAngle = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
            rcs(i)= rcsellipsoid(thickness,thickness,B_sum_sqrt/2,0,ThetaAngle);
        end
        amp=sqrt(rcs);
        returned_signal=returned_signal+amp.*Signal(time-2*distance/c);
    end
end

%% calculate the signal returned from the wheels
for j=1:2 % front and back wheels
    for k=1:2 % left and right wheels
        message=['calculating the (',int2str(j),',',int2str(k),') wheel'];
        disp(message);
        for m=1:rings
            for n=1:points
                P1=wheelPositionData(j,k).points(m,n).location;
                if n==points
                    P2=wheelPositionData(j,k).points(m,1).location;
                else
                    P2=wheelPositionData(j,k).points(m,n+1).location;
                end
                P=(P1+P2)/2;
                dx=P(1,:)-radar_x;
                dy=P(2,:)-radar_y;
                dz=P(3,:)-radar_z;
                distance=sqrt(dx.^2+dy.^2+dz.^2);
                % calculate the RCS
                rcs=zeros(size(time));
                for i=1:length(time)
                    A=[dx(i),dy(i),dz(i)];
                    B=[P1(1,i)-P2(1,i),P1(2,i)-P2(2,i),P1(3,i)-P2(3,i)];
                    A_dot_B = dot(A,B);
                    A_sum_sqrt = sqrt(sum(A.*A));
                    B_sum_sqrt = sqrt(sum(B.*B));
                    ThetaAngle = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
                    rcs(i)= rcsellipsoid(thickness,thickness,B_sum_sqrt/2,0,ThetaAngle);
                end
                amp=sqrt(rcs);
                returned_signal=returned_signal+amp.*Signal(time-2*distance/c);
            end
        end
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

