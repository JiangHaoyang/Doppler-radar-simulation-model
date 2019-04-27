function [ output ] = RectangleWave( time, width )
%generate a rectangle wave, whose size is the same as 'time'. When time >=
%0 and <= width, output = 1. Otherwise, 0;
output=zeros(size(time));
output(time>=0 & time<=width)=1;
end

