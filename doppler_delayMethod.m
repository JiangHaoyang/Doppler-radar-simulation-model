function [ DopplerDelayAnimation,spectrogram,spectrogram_timeAxis,spectrogram_dopplerAxis,spectrogramCmplx ] = doppler_delayMethod( Signal,shortTimeParameters,longTimeParameters,...
    carParameters,wheelPositionData,radarLocation,controlFlag)
% Signal is a function: Signal=@(t) ...
% shortTimeParameters=[sample_rate,startF,endF,Signal_duration,Pulse_duration,c];
% longTimeParameters=[startTime,endTime,sampleRate];
% carParameters=[carLength,carWidth,radius,thickness,ringsInWheel,pointsInRing];
% controlFlag=[generateRetuenSignal,generateFrames];

warning off backtrace;
generateRetuenSignal=controlFlag(1);
generateFrames=controlFlag(2);
%% set parameters
startTime=longTimeParameters(1);
endTime=longTimeParameters(2);
longSampleRate=longTimeParameters(3);
dT=1/longSampleRate;
time=startTime:dT:endTime;

shortSample_rate=shortTimeParameters(1);
startF=shortTimeParameters(2);
endF=shortTimeParameters(3);
Signal_duration=shortTimeParameters(4);
Pulse_duration=shortTimeParameters(5);
c=shortTimeParameters(6);
Pulse_time=0:1/shortSample_rate:Signal_duration;


thickness=carParameters(4);
rings=carParameters(5);
points=carParameters(6);

radar_x=radarLocation(1);
radar_y=radarLocation(2);
radar_z=radarLocation(3);

current_filepath = mfilename('fullpath');
current_filepath=fileparts(current_filepath);
exampleParameter_path=[current_filepath,'\Sample_1\Parameters'];
mkdir(exampleParameter_path);
save([exampleParameter_path,'\','parameters'],'longSampleRate','shortTimeParameters');

%% generate raw data
% set the data save path
if generateRetuenSignal=='y'
    current_filepath = mfilename('fullpath');
    current_filepath=fileparts(current_filepath);
    retuened_signal_data_path=[current_filepath,'\Sample_1\ReturnSignal'];
    if isdir(retuened_signal_data_path)
        warning('pulses data already exist in the direction');
        userInput=input('Do you want to clear the old pulses data and continue?(yes/no)\n','s');
        if strcmp(userInput,'yes')
            warning('previous pulses data are cleared');
            rmdir(retuened_signal_data_path,'s');
        else
            warning('system stop');
            generateFrames='n';
            generateRetuenSignal='n';
        end
    end
end

if generateRetuenSignal=='y'
    % make file direction
    mkdir(retuened_signal_data_path);
    % calculate the pulses one by one and save
    
    % parameters for displaying process
    old_workingRate=-1;
    new_workingRate=0;
    
    for i=1:length(time)
        returned_pulse=zeros(size(Pulse_time));
        for j=1:2
            for k=1:2
                % signals from axles
                AxleIndex=(j-1)*2+k;
                P1=wheelPositionData(AxleIndex).center(:,i);
                if AxleIndex==4
                    P2=wheelPositionData(1).center(:,i);
                else
                    P2=wheelPositionData(AxleIndex+1).center(:,i);
                end
                P=(P1+P2)/2;
                dx=P(1)-radar_x;
                dy=P(2)-radar_y;
                dz=P(3)-radar_z;
                distance=sqrt(dx.^2+dy.^2+dz.^2);
                % calculate the RCS
                A=[dx,dy,dz];
                B=[P1(1)-P2(1),P1(2)-P2(2),P1(3)-P2(3)];
                A_dot_B = dot(A,B);
                A_sum_sqrt = sqrt(sum(A.*A));
                B_sum_sqrt = sqrt(sum(B.*B));
                ThetaAngle = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
                rcs = rcsellipsoid(thickness,thickness,B_sum_sqrt/2,0,ThetaAngle);
                amp = sqrt(rcs);
                % returned signal
                % S_r(t)=S(t-2*distance/c)
                returned_pulse=returned_pulse+amp*Signal(Pulse_time-2*distance/c);
                
                % signals from wheels
                for m=1:rings
                    for n=1:points
                        P1=wheelPositionData(j,k).points(m,n).location(:,i);
                        if n==points
                            P2=wheelPositionData(j,k).points(m,1).location(:,i);
                        else
                            P2=wheelPositionData(j,k).points(m,n+1).location(:,i);
                        end
                        P=(P1+P2)/2;
                        dx=P(1)-radar_x;
                        dy=P(2)-radar_y;
                        dz=P(3)-radar_z;
                        distance=sqrt(dx.^2+dy.^2+dz.^2);
                        % calculate the RCS
                        A=[dx,dy,dz];
                        B=[P1(1)-P2(1),P1(2)-P2(2),P1(3)-P2(3)];
                        A_dot_B = dot(A,B);
                        A_sum_sqrt = sqrt(sum(A.*A));
                        B_sum_sqrt = sqrt(sum(B.*B));
                        ThetaAngle = acos(A_dot_B ./ (A_sum_sqrt .* B_sum_sqrt));
                        rcs = rcsellipsoid(thickness,thickness,B_sum_sqrt/2,0,ThetaAngle);
                        amp = sqrt(rcs);
                        % returned signal
                        % S_r(t)=S(t-2*distance/c)
                        returned_pulse=returned_pulse+amp*Signal(Pulse_time-2*distance/c);
                    end
                end
            end
        end
        file_name = sprintf('Pulses%08i.mat',i);
        save([retuened_signal_data_path,'\',file_name],'returned_pulse');
        
        % Display the process
        step=2;
        new_workingRate=floor(i/length(time)*100/step);
        if(new_workingRate>old_workingRate)
            fprintf('Step 1/2, generating raw pulses data(%i%%)\n',new_workingRate*step);
            old_workingRate=new_workingRate;
        end
    end
end
%% generate frames data
if generateFrames=='y';
    
    % set the data save direction
    current_filepath = mfilename('fullpath');
    current_filepath=fileparts(current_filepath);
    frames_RawData_data_path=[current_filepath,'\Sample_1\FrameData'];
    retuened_signal_data_path=[current_filepath,'\Sample_1\ReturnSignal'];
    if isdir(frames_RawData_data_path)
        warning('frames data already exist in the direction');
        userInput=input('Do you want to clear the old frames data and continue?(yes/no)\n','s');
        if strcmp(userInput,'yes')
            warning('previous frames data are cleared');
            rmdir(frames_RawData_data_path,'s');
        else
            generateFrames='n';
            warning('system stop');
        end
    end
end

if generateFrames=='y';
    % make file direction
    mkdir(frames_RawData_data_path);
    
    plusesNumber=100;% numbers of pulses per frame
    increment=16;% increment of next frame
    pulses_num=length(time);
    frames_num=length(plusesNumber:increment:pulses_num);
    
    % impulseRespon
    ImpResp = Signal(Pulse_time);
    
    FrameData=zeros(plusesNumber,2*length(ImpResp)-1);
    
    % doppler_point
    doppler_point=plusesNumber;
    m = nextpow2(4*doppler_point);
    doppler_point=2^m;
    
    
    % save axis data;
    delay_axis=linspace(-Signal_duration,Signal_duration,2*length(ImpResp)-1);
    doppler_axis=((1:doppler_point)-1)./doppler_point*longSampleRate;
    doppler_axis=doppler_axis-doppler_axis(ceil((length(doppler_axis)+1)/2));
    distance_axis=c*delay_axis/2;
    save([frames_RawData_data_path,'\','doppler_delay_image_axis.mat'],'delay_axis','doppler_axis','distance_axis');
    
    % parameters for displaying process
    old_workingRate=-1;
    new_workingRate=0;
    
    for i=1:(plusesNumber-increment)
        file_name = sprintf('Pulses%08i.mat',i);
        load([retuened_signal_data_path,'\',file_name]);
        FrameData(increment+i,:)=xcorr(returned_pulse,ImpResp);
        
        % Display the process
        step=5;
        new_workingRate=floor(i/(plusesNumber-increment)*100/step);
        if(new_workingRate>old_workingRate)
            fprintf('Step 2/2, loading data(%i%%)\n',new_workingRate*step);
            old_workingRate=new_workingRate;
        end
    end
    
    % empty animation file
    DopplerDelayAnimation(frames_num) = struct('cdata',[],'colormap',[]);
    % empty spectrum image
    spectrogram=zeros(doppler_point,frames_num);
    spectrogramCmplx=zeros(doppler_point,frames_num);
    
    % parameters for displaying process
    old_workingRate=-1;
    new_workingRate=0;
    
    myfig=figure('name','doppler_delay image');
    for i=1:frames_num
        
        FrameData(1:(plusesNumber-increment),:)=FrameData(increment+1:plusesNumber,:);
        for m=1:increment
            file_name = sprintf('Pulses%08i.mat',(plusesNumber-2*increment)+i*increment+m);
            load([retuened_signal_data_path,'\',file_name],'returned_pulse');
            FrameData((plusesNumber-increment)+m,:)=xcorr(returned_pulse,ImpResp);
        end
        RawData=fft(FrameData,doppler_point,1);
        file_name = sprintf('Frame%04i.mat',i);
        save([frames_RawData_data_path,'\',file_name],'RawData','-v7.3');
        clf(myfig);
        % draw and save animation frames
        Image_data=fftshift(20*log10(abs(RawData)),1);
        max_v=max(max(Image_data));
        imagesc(delay_axis,doppler_axis,Image_data,[max_v-50,max_v]);
        colormap(jet);
        axis([delay_axis(1),delay_axis(end),doppler_axis(1),doppler_axis(end)]);
        xlabel('delay');
        ylabel('doppler');
        DopplerDelayAnimation(i)=getframe(myfig);
        drawnow;
        
        % generate the spectrum
        spectrogram(:,i)=sum(abs(RawData),2);
        spectrogramCmplx(:,i)=sum(RawData,2);
        % Display the process
        step=1;
        new_workingRate=floor(i/frames_num*100/step);
        if(new_workingRate>old_workingRate)
            fprintf('Step 2/2, processing frames data(%i%%)\n',new_workingRate*step);
            old_workingRate=new_workingRate;
        end
        
    end
    
    % save spectrum axis data;
    spectrogram_timeAxis=time(plusesNumber:increment:pulses_num);
    spectrogram_dopplerAxis=doppler_axis;
    save([frames_RawData_data_path,'\','spectrogram_axis'],'spectrogram_timeAxis','spectrogram_dopplerAxis');
    save([frames_RawData_data_path,'\','Frames_Animation'],'DopplerDelayAnimation');
    save([frames_RawData_data_path,'\','Specturogram'],'spectrogram','spectrogramCmplx');
else
    % assign null value to all returned parameters
    DopplerDelayAnimation=0;
    spectrogram=0;
    spectrogram_timeAxis=0;
    spectrogram_dopplerAxis=0;
    spectrogramCmplx=0;
end
warning on backtrace;
end

