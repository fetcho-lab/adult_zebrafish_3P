%makes an avi movie from a red/green tif image series
clear all, close all
framerate=18;


[tifFile, path] = uigetfile('*.tif','Select .tif file to convert to movie');
cd(path);

imgdata = readImage(tifFile);
img_min = 650;
img_max = 1596;

imgdata = [double(imgdata) - img_min]/(img_max - img_min);
imgdata(imgdata<0) = 0;
imgdata(imgdata>1) = 1;
imgdata = single(imgdata);


% outVideo = VideoWriter([caFile(1:end-4),'.avi'],'Motion JPEG AVI');
outVideo = VideoWriter([tifFile(1:end-4),'.avi'],'Uncompressed AVI');
set(outVideo,'FrameRate',framerate);
open(outVideo);


for t=1:size(imgdata,3)
    
    
    writeVideo(outVideo,imgdata(:,:,t));
    if mod(outVideo.FrameCount,100) == 0
        disp(sprintf('%3.0f frames written',outVideo.FrameCount));
    end
end
close(outVideo);
