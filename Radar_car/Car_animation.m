function [ animation ] = Car_animation( wheelPositionData,carParameters,radarLocation,AXIS )
%carParameters=[carLength,carWidth,radius,thickness,ringsInWheel,pointsInRing];

% get the time length
time=size(wheelPositionData(1,1).points(1,1).location);
time=time(2);

thickness=carParameters(4);
rings=carParameters(5);
points=carParameters(6);

playSpeed=200;
index=0;
myfig=figure('name','animation');
for i=1:1*playSpeed:time
    clf(myfig);
    hold on;
    % draw the Radar
    [X,Y,Z]=ellipsoid(radarLocation(1),radarLocation(2),radarLocation(3),0.5,0.5,0.5,5);
    surf(X,Y,Z);
    % draw the fornt axle
    P1=wheelPositionData(1,1).center(:,i);P1=P1';
    P2=wheelPositionData(1,2).center(:,i);P2=P2';
    semi_radius=0.5*sqrt(sum((P1-P2).^2));
    [X,Y,Z] = ellipsoid2P(P1,P2,thickness,thickness,semi_radius,10);
    surf(X,Y,Z);
    
    % draw the back axle
    P1=wheelPositionData(2,1).center(:,i);P1=P1';
    P2=wheelPositionData(2,2).center(:,i);P2=P2';
    semi_radius=0.5*sqrt(sum((P1-P2).^2));
    [X,Y,Z] = ellipsoid2P(P1,P2,thickness,thickness,semi_radius,10);
    surf(X,Y,Z);
    
    % draw the left axle
    P1=wheelPositionData(1,1).center(:,i);P1=P1';
    P2=wheelPositionData(2,1).center(:,i);P2=P2';
    semi_radius=0.5*sqrt(sum((P1-P2).^2));
    [X,Y,Z] = ellipsoid2P(P1,P2,thickness,thickness,semi_radius,10);
    surf(X,Y,Z);
    
    % draw the right axle
    P1=wheelPositionData(1,2).center(:,i);P1=P1';
    P2=wheelPositionData(2,2).center(:,i);P2=P2';
    semi_radius=0.5*sqrt(sum((P1-P2).^2));
    [X,Y,Z] = ellipsoid2P(P1,P2,thickness,thickness,semi_radius,10);
    surf(X,Y,Z);
    
    % draw wheels
    for j=1:2
        for k=1:2
            for m=1:rings
                for n=1:points
                    P1=wheelPositionData(j,k).points(m,n).location(:,i);P1=P1';
                    if n+1>points
                        P2=wheelPositionData(j,k).points(m,1).location(:,i);P2=P2';
                    else
                        P2=wheelPositionData(j,k).points(m,n+1).location(:,i);P2=P2';
                    end
                    semi_radius=0.5*sqrt(sum((P1-P2).^2));
                    [X,Y,Z] = ellipsoid2P(P1,P2,thickness,thickness,semi_radius,10);
                    surf(X,Y,Z);
                end
            end
        end
    end
    index=index+1;
    xlabel('x');ylabel('y');zlabel('z');
    axis equal;
    axis(AXIS);
    view(70,10);
    drawnow;
    animation(index)=getframe(myfig);
end
end

