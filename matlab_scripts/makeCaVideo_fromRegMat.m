%makes an avi movie from a red/green tif image series
clear all, close all
framerate=4.24*12;



[caFile,path] = uigetfile('*.mat','Select .mat file to convert to movie');
cd(path)
load(caFile)

imgRange = double( [min(imgdata(:)), max(imgdata(:))] );

% outVideo = VideoWriter([caFile(1:end-4),'.avi'],'Motion JPEG AVI');
outVideo = VideoWriter([caFile(1:end-4),'.avi'],'Uncompressed AVI');
open(outVideo);

imgMin = 0;
imgMax = imgRange(2)/2;

window_halfsize = 3;

%compress time series

% for t=1:2:size(imgdata,3)
%     imgdata_compressed(:,:,(t+1)/2) = mean(imgdata(:,:,t:t+1),3);
% end
% imgdata = imgdata_compressed;

for t=1:size(imgdata,3)
    negCon = max([1, t-window_halfsize]);
    posCon = min([size(imgdata,3), t+window_halfsize]);
    
    cFrame = mean( imgdata(:,:,negCon:posCon) , 3);
%     cFrame = medfilt2(cFrame);
    
    cFrame = (cFrame-imgMin)./imgMax;
    cFrame(cFrame<0) = 0;
    cFrame = uint8(cFrame*255);
    
    
    writeVideo(outVideo,cFrame);
    if mod(outVideo.FrameCount,100) == 0
        disp(sprintf('%3.0f frames written',outVideo.FrameCount));
    end
end
close(outVideo);
