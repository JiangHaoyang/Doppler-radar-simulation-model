function [ DopplerDelayAnimation,spectrogram,spectrogram_timeAxis,spectrogram_dopplerAxis,spectrogramCmplx ] = ...
    Doppler_delayMethod_forWindTurbine( Signal,shortTimeParameters,longTimeParameters,...
    model_data,radarLocation,controlFlag)
% Signal is a function: Signal=@(t) ...
% shortTimeParameters=[sample_rate,startF,endF,Signal_duration,Pulse_duration,c];
% longTimeParameters=[startTime,endTime,sampleRate];
% controlFlag=[generateRetuenSignal,generateFrames];
%
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
%% function execution
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

radar_x=radarLocation(1);
radar_y=radarLocation(2);
radar_z=radarLocation(3);

current_filepath = mfilename('fullpath');
current_filepath=fileparts(current_filepath);
exampleParameter_path=[current_filepath,'\Sample_1\Parameters'];
mkdir(exampleParameter_path);
modelparameter.longSampleRate=longSampleRate;
modelparameter.shortTimeParameters=shortTimeParameters;
save([exampleParameter_path,'\','parameters'],'modelparameter');

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
    
    % generate fixed returned signal, such as signal returned from the tower,
    % as the radar position is fix and the tower can not move, the returned
    % signal are always the same
    TMP=zeros(size(Pulse_time));
    % signal from tower
    for k=1:3
        P1=model_data.tower(1,:,k);
        P2=model_data.tower(2,:,k);
        P=(P1+P2)/2;
        dx=P(1)-radar_x;
        dy=P(2)-radar_y;
        dz=P(3)-radar_z;
        distance=sqrt(dx.^2+dy.^2+dz.^2);
        rcs=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% not write
        amp=sqrt(rcs);
        TMP=TMP+amp.*Signal(Pulse_time-2*distance/c);
    end
    % signal from nacelle
    P=model_data.nacelle(1,:);
    dx=P(1)-radar_x;
    dy=P(2)-radar_y;
    dz=P(3)-radar_z;
    distance=sqrt(dx.^2+dy.^2+dz.^2);
    rcs=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% not write
    amp=sqrt(rcs);
    TMP=TMP+amp.*Signal(Pulse_time-2*distance/c);
    % signal from roter
    P=model_data.rotor(1,:);
    dx=P(1)-radar_x;
    dy=P(2)-radar_y;
    dz=P(3)-radar_z;
    distance=sqrt(dx.^2+dy.^2+dz.^2);
    rcs=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% not write
    amp=sqrt(rcs);
    TMP=TMP+amp.*Signal(Pulse_time-2*distance/c);
    
    % parameters for displaying process
    old_workingRate=-1;
    new_workingRate=0;
    % signal from blades
    for i=1:length(time)
        returned_pulse=TMP;
        for m=1:3
            for n=1:4
                P1=model_data.bladePosition(:,i,m,n);
                P2=model_data.bladePosition(:,i,m,n+1);
                P=(P1+P2)/2;
                dx=P(1,:)-radar_x;
                dy=P(2,:)-radar_y;
                dz=P(3,:)-radar_z;
                distance=sqrt(dx.^2+dy.^2+dz.^2);
                rcs=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% not write
                amp=sqrt(rcs);
                returned_pulse=returned_pulse+amp.*Signal(Pulse_time-2*distance/c);
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
    
    plusesNumber=200;% numbers of pulses per frame
    increment=50;% increment of next frame
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

