function SaveAVI( animation )
[~,Frames_num]=size(animation);
writerObj=VideoWriter('movie.avi');
writerObj.FrameRate=20;
open(writerObj);
myfig=figure('name','animation');
for i=1:Frames_num
    disp(i);
    clf(myfig);
    imshow(animation(i).cdata);
    F=getframe(myfig);
    writeVideo(writerObj,F);
end
close(writerObj);
end
