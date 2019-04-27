function [ model_data ] = WindTurbineModel( ModelParameters,timeParameters )
% input parameters:
% timeParameters=[startTime,endTime,sampleRate];
% ModelParameters=[tower1_height , tower1_radius  , 0 , 0 ;        % bottom of the tower
%                  tower2_height , tower2_radius  , 0 , 0 ;        % middle of the tower
%                  tower3_height , tower3_radius  , 0 , 0 ;        % top of the tower
%                  nacelle_x     , nacelle_y,     ,nacelle_z ,0 ;  % size of the nacelle
%                  nacelle_bias_x, nacelle_bias_y , 0 , 0 ;        % tower Top At The Bottom Of TheNacelle
%                  rotor_radius  , rotorCneterBeyondTheNacelle,0 ,0 ;
%                  blade_1_height, blade_1_radius ,0 , 0 ;
%                  blade_2_x     , blade_2_y      ,blade_2_z , blade_2_rotate;
%                  blade_3_x     , blade_3_y      ,blade_3_z , blade_3_rotate;
%                  blade_4_x     , blade_4_y      ,blade_4_z , blade_4_rotate;
%                  RotateSpeed   , 0 , 0 , 0 ];

% output parameters:
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
%% extract information from the WindTurbineModel
tower1_height=ModelParameters(1,1);
tower1_radius=ModelParameters(1,2);
tower2_height=ModelParameters(2,1);
tower2_radius=ModelParameters(2,2);
tower3_height=ModelParameters(3,1);
tower3_radius=ModelParameters(3,2);
nacelle_x=ModelParameters(4,1);
nacelle_y=ModelParameters(4,2);
nacelle_z=ModelParameters(4,3);
nacelle_bias_x=ModelParameters(5,1);
% nacelle_bias_y=ModelParameters(5,2);
rotor_radius=ModelParameters(6,1);
rotorCneterBeyondTheNacelle=ModelParameters(6,2);
bladeInfo=ModelParameters(7:10,:);% bladeInfo contains all information of the shape of the blades
RotateSpeed=ModelParameters(11,1);
%% set time
startTime=timeParameters(1);
endTime=timeParameters(2);
sampleRate=timeParameters(3);
dT=1/sampleRate;
time=startTime:dT:endTime;

%% calculate coordinate
% tower
tower_H1=[0,0,0];
tower_H2=[0,0,tower1_height];
tower_H3=[0,0,tower1_height+tower2_height];
tower_H4=[0,0,tower1_height+tower2_height+tower3_height];
model_data.tower=zeros(3,3,3);
model_data.tower(:,:,1)=[tower_H1;tower_H2;tower1_height,tower1_radius,0];
model_data.tower(:,:,2)=[tower_H2;tower_H3;tower2_height,tower2_radius,0];
model_data.tower(:,:,3)=[tower_H3;tower_H4;tower3_height,tower3_radius,0];

% nacelle
nacelleCenter=[nacelle_bias_x-nacelle_x/2,0,tower_H4(3)+nacelle_z/2];
model_data.nacelle=[nacelleCenter;nacelle_x,nacelle_y,nacelle_z];

% rotor
rotorCenter=[(nacelleCenter(1)+nacelle_x/2+rotorCneterBeyondTheNacelle),0,tower_H4(3)+nacelle_z/2];
model_data.rotor=[rotorCenter;rotor_radius,0,0];

% blades radius
bladeRadius=zeros(1,5);
bladeRadius(1)=rotor_radius;
bladeRadius(2)=bladeRadius(1)+bladeInfo(1,1);
bladeRadius(3)=bladeRadius(2)+bladeInfo(2,1);
bladeRadius(4)=bladeRadius(3)+bladeInfo(3,1);
bladeRadius(5)=bladeRadius(4)+bladeInfo(4,1);
model_data.bladeInfo=bladeInfo;
% blades angle
iniAngle=zeros(1,3);
iniAngle(1)=0;%0 degree
iniAngle(2)=2*pi/3;% 120 degree
iniAngle(3)=2*pi*2/3;% 240 degree
% blades
model_data.bladePosition=zeros(3,length(time),3,5);
for i=1:length(time)
    rotorRotateAngle=RotateSpeed*2*pi/60*time(i);
    for m=1:3 % 3 blades
        Angle=rotorRotateAngle+iniAngle(m);
        % if sin(Angle)=0 the system will crashes
        if(sin(Angle)==0)
            Angle=Angle+0.000001;
        end
        for n=1:5 % P1 to P5
            model_data.bladePosition(:,i,m,n)=bladeRadius(n)*[0;sin(Angle);cos(Angle)]+rotorCenter';
        end
    end
end

end

















