zdclear; close all;
% a script to extract ca/functional imaging data from 3P fluorescence
% channel, and generate registered .tif stack
[fSFile, fSPath] = uigetfile('*.tif', 'Select Image File');

cd(fSPath);

fDat = imfinfo(fSFile);
img0 = [];

offsets = [];
nChannelsAcq = 1;
%average group size  ==  5 or 10. then use to generate an interpolation
%line
GroupSize = 10;

imagineSeries = zeros(fDat(1).Height,fDat(1).Width,numel(fDat));
indices = [1:nChannelsAcq:numel(fDat)];

parfor k = 1:numel(fDat)/nChannelsAcq
    index = indices(k);
    dataImage =  double( imread(fSFile,index) );
    dataImage = dataImage - median(dataImage(:));
    dataImage(dataImage<0) = 0;
    
    imagineSeries(:,:,k) = dataImage;
end

for k=1:size(imagineSeries,3)/GroupSize
    groupImagination(:,:,k) = mean(imagineSeries(:,:,(k-1)*GroupSize+1:k*GroupSize),3);
end

offsets = [];
img0 = groupImagination(:,:,1);

parfor k = 2:size(groupImagination,3)
%read even numbers only for fluorescence channel

img = groupImagination(:,:,k);
% img = medfilt2(img,'symmetric');

imgX = normxcorr2(img0,img);

% offset found by correlation
maxCorr = max(imgX(:));
[aY,aX] = find(imgX == maxCorr);

offsets(k-1,:) = [aX,aY]-fDat(1).Height;

end

timeAx = GroupSize:GroupSize:numel(fDat);
goodOffsets = sum(offsets,2) < 10;

%calculate a drift rate based on this. 
driftRate_aX = timeAx(goodOffsets)'\offsets(goodOffsets,2);
% driftRate_aY= timeAx(1:end-1)'\offsets(:,1);
driftRate_aY= timeAx(goodOffsets)'\offsets(goodOffsets,1);

imgCorrectionaX = -[1:numel(fDat)]*driftRate_aX;
imgCorrectionaY = -[1:numel(fDat)]* driftRate_aY;

parfor k= 1:numel(fDat)
    X = imgCorrectionaX(k);
    Y = imgCorrectionaY(k);
    imgCorrected(:,:,k) = imtranslate(imagineSeries(:,:,k), [Y,X]);
    
end

figure(1),
% subplot(3,2,1)
% imagesc(groupImagination(:,:,1));
% title('Raw Frame 1');
% 
% subplot(3,2,2)
% imagesc(imgCorrected(:,:,1));
% title('Motion Corrected Frame 1');
% 
% subplot(3,2,3)
% imagesc(groupImagination(:,:,end));
% title('Raw Last Frame');
subplot(3,2,1)
imagesc(groupImagination(:,:,end));
title('Last frame uncorrected');

subplot(3,2,2);
imagesc(mean(imgCorrected(:,:,end-GroupSize:end),3));
title('Motion Corrected Last Frame');

subplot(3,2,3);
imagesc(mean(groupImagination,3));
title('Average of Raw');

subplot(3,2,4);
imagesc(mean(imgCorrected,3));
title('Average of Corrected');

subplot(3,2,5);
imagesc(groupImagination(:,:,end) - groupImagination(:,:,1));
title('Difference Raw');

subplot(3,2,6)
imagesc(mean(imgCorrected(:,:,end-GroupSize:end),3) - mean(imgCorrected(:,:,1:GroupSize),3));
title('Difference Motion Corrected');

disp('Press enter to save');
pause;

mkdir('Registered');
fs = filesep;
savefig(['Registered',fs,fSFile(1:end-4),'.fig']);
imgdata = uint16(imgCorrected );
save(['Registered',fs,fSFile(1:end-4),'.mat'],'imgdata','-v7.3');

for frame=1:size(imgdata, 3)

    tiffObj = Tiff(['Registered',fs,fSFile(1:end-4),'_driftreg.tif'],'a');

    tiffObj.setTag('Photometric',Tiff.Photometric.LinearRaw);
    tiffObj.setTag('BitsPerSample',16);
    tiffObj.setTag('ImageWidth',size(imgdata,2));
    tiffObj.setTag('ImageLength',size(imgdata,1));
    tiffObj.setTag('SamplesPerPixel',1);
    tiffObj.setTag('Compression',Tiff.Compression.PackBits);
    tiffObj.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);

    tiffObj.write( imgdata(:,:,frame) );        


    if frame == size(imgdata, 3)
        writeDirectory(tiffObj);
        tiffObj.close();
    end
end
