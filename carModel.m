function wheelPositionData = carModel( carParameters,timeParameters,trajectory )
%carParameters=[carLength,carWidth,radius,thickness,ringsInWheel,pointsInRing];
%timeParameters=[startTime,endTime,sampleRate];
%trajectory is a funciton

carLength=carParameters(1);
carWidth=carParameters(2);
radius=carParameters(3);
% thickness=carParameters(4);
rings=carParameters(5);
points=carParameters(6);

startTime=timeParameters(1);
endTime=timeParameters(2);
sampleRate=timeParameters(3);
dT=1/sampleRate;
time=startTime:dT:endTime;

% car center coordinate
car_center=zeros(3,length(time));
car_center(1,:)=trajectory(time);

% frontLeft wheel center
temp=zeros(3,length(time));
temp(1,:)=carLength/2;
temp(2,:)=carWidth/2;
FLWheel_center=car_center+temp;

% frontRignt wheel center
temp=zeros(3,length(time));
temp(1,:)=carLength/2;
temp(2,:)=-carWidth/2;
FRWheel_center=car_center+temp;

% backLeft wheel center
temp=zeros(3,length(time));
temp(1,:)=-carLength/2;
temp(2,:)=carWidth/2;
BLWheel_center=car_center+temp;

% backRight wheel center
temp=zeros(3,length(time));
temp(1,:)=-carLength/2;
temp(2,:)=-carWidth/2;
BRWheel_center=car_center+temp;


travel_distance=car_center(1,:)-car_center(1,1);
% data for front left wheel
FLWheel(rings,points).location=zeros(3,length(time));
for m=1:rings
    m_radius=m/rings*radius;
    for n=1:points
        ini_angel=2*pi*n/points;
        location=zeros(3,length(time));
        location(1,:)=m_radius*cos(travel_distance/radius+ini_angel);
        location(3,:)=-m_radius*sin(travel_distance/radius+ini_angel);
        FLWheel(m,n).location=location+FLWheel_center;
    end
end

% data for front right wheel
FRWheel(rings,points).location=zeros(3,length(time));
for m=1:rings
    m_radius=m/rings*radius;
    for n=1:points
        ini_angel=2*pi*n/points;
        location=zeros(3,length(time));
        location(1,:)=m_radius*cos(travel_distance/radius+ini_angel);
        location(3,:)=-m_radius*sin(travel_distance/radius+ini_angel);
        FRWheel(m,n).location=location+FRWheel_center;
    end
end

% data for back left wheel
BLWheel(rings,points).location=zeros(3,length(time));
for m=1:rings
    m_radius=m/rings*radius;
    for n=1:points
        ini_angel=2*pi*n/points;
        location=zeros(3,length(time));
        location(1,:)=m_radius*cos(travel_distance/radius+ini_angel);
        location(3,:)=-m_radius*sin(travel_distance/radius+ini_angel);
        BLWheel(m,n).location=location+BLWheel_center;
    end
end

% data for back right wheel
BRWheel(rings,points).location=zeros(3,length(time));
for m=1:rings
    m_radius=m/rings*radius;
    for n=1:points
        ini_angel=2*pi*n/points;
        location=zeros(3,length(time));
        location(1,:)=m_radius*cos(travel_distance/radius+ini_angel);
        location(3,:)=-m_radius*sin(travel_distance/radius+ini_angel);
        BRWheel(m,n).location=location+BRWheel_center;
    end
end

wheelPositionData(1,1).name='front_left';
wheelPositionData(1,1).center=FLWheel_center;
wheelPositionData(1,1).points=FLWheel;
wheelPositionData(1,2).name='front_right';
wheelPositionData(1,2).center=FRWheel_center;
wheelPositionData(1,2).points=FRWheel;
wheelPositionData(2,1).name='back_left';
wheelPositionData(2,1).center=BLWheel_center;
wheelPositionData(2,1).points=BLWheel;
wheelPositionData(2,2).name='back_right';
wheelPositionData(2,2).center=BRWheel_center;
wheelPositionData(2,2).points=BRWheel;

end

