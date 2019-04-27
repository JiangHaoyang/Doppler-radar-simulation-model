function [ animation ] = WindTurbine_animation( windTurbine_modelData,radarLocation,AXIS )
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
time=size(windTurbine_modelData.bladePosition);
time=time(2);
playSpeed=50;
index=0;
myfig=figure('name','animation');
for i=1:playSpeed:time
    clf(myfig);hold on;
    % draw the tower
    for k=1:3
        P1=windTurbine_modelData.tower(1,:,k);
        P2=windTurbine_modelData.tower(2,:,k);
        h=windTurbine_modelData.tower(3,1,k);
        r=windTurbine_modelData.tower(3,2,k);
        [X,Y,Z]=cylinder2P(P1,P2,h,r,20);
        surf(X,Y,Z);
    end
    
    % draw the nacelle
    nacelleCenter=windTurbine_modelData.nacelle(1,:);
    nacelleSize=windTurbine_modelData.nacelle(2,:);
    [X,Y,Z]=cube(nacelleCenter,nacelleSize(1),nacelleSize(2),nacelleSize(3));
    surf(X,Y,Z);
    
    % draw the rotor
    rotorCenter=windTurbine_modelData.rotor(1,:);
    rotorRadius=windTurbine_modelData.rotor(2,1);
    P1=rotorCenter;P1(3)=P1(3)+rotorRadius;
    P2=rotorCenter;P2(3)=P2(3)-rotorRadius;
    [X,Y,Z]=ellipsoid2P(P1,P2,rotorRadius,rotorRadius,rotorRadius,10);
    surf(X,Y,Z);
    
    % draw the blades
    for m=1:3 % 3 blades
        for n=1:4 % 4 parts of each blade
            P1=windTurbine_modelData.bladePosition(:,i,m,n);
            P2=windTurbine_modelData.bladePosition(:,i,m,n+1);
            if (n==1)
                h=windTurbine_modelData.bladeInfo(1,1);
                r=windTurbine_modelData.bladeInfo(1,2);
                [X,Y,Z]=cylinder2P(P1',P2',h,r,20);
            else
                ellipsoid_x=windTurbine_modelData.bladeInfo(n,1);
                ellipsoid_y=windTurbine_modelData.bladeInfo(n,2);
                ellipsoid_z=windTurbine_modelData.bladeInfo(n,3);
                ellipsoid_angle=windTurbine_modelData.bladeInfo(n,4);
                [X,Y,Z]=ellipsoid2P1A(P1',P2',ellipsoid_angle,ellipsoid_x/2,ellipsoid_y/2,ellipsoid_z/2,20);
            end
            surf(X,Y,Z);
        end
    end
    % draw the radar
    radarRadius=3;
    P1=radarLocation;P1(3)=P1(3)+radarRadius;
    P2=radarLocation;P2(3)=P2(3)-radarRadius;
    [X,Y,Z]=ellipsoid2P(P1,P2,radarRadius,radarRadius,radarRadius,10);
    surf(X,Y,Z);
    % draw the targeting vector
    quiver3(radarLocation(1),radarLocation(2),radarLocation(3),rotorCenter(1)-radarLocation(1),rotorCenter(2)-radarLocation(2),rotorCenter(3)-radarLocation(3),'k','filled','LineWidth',1);
    xlabel('x(m)');
    ylabel('y(m)');
    zlabel('z(m)');
    axis equal;
    axis(AXIS);
    view(140,20);
    drawnow;
    
    % save the current fig
    index=index+1;
    animation(index)=getframe(myfig);
end
end


















