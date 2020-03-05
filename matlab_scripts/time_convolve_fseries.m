clear; close all;
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
halfGroupSize = 3;

imagineSeries = zeros(fDat(1).Height,fDat(1).Width,numel(fDat));
indices = [1:nChannelsAcq:numel(fDat)];

parfor k = 1:numel(fDat)/nChannelsAcq
    index = indices(k);
    dataImage =  double( imread(fSFile,index) );
    dataImage = dataImage - median(dataImage(:));
    dataImage(dataImage<0) = 0;
    
    imagineSeries(:,:,k) = dataImage;
end

for k=1:size(imagineSeries,3)
    time0 = max([1, k-halfGroupSize]);
    time1 = min([size(imagineSeries,3), k+halfGroupSize]);
    groupImagination(:,:,k) = mean(imagineSeries(:,:,time0:time1),3);
end


mkdir('Convolved');
fs = filesep;
imgdata = uint16(groupImagination );
%save(['Convolved',fs,fSFile(1:end-4),'_convolved.mat'],'imgdata','-v7.3');

for frame=1:size(imgdata, 3)

    tiffObj = Tiff(['Convolved',fs,fSFile(1:end-4),'_convolved.tif'],'a');

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
