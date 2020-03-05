clear; 
tiffPath = uigetdir(pwd,'Select Source Data for Related Stacks');
cd(tiffPath);
%for this to work: put all related stacks in the same directory. They must
%have consistent motor coordinates from the software (I think) for this to
%work right, and must have the same step size in Z. 

%+Y, fov left (that is image is moving to the right)
%+X, fov down (due to mirror), so image is moving up

fullFieldMicrons = 400; %this value is 400 on Dimitre's rig and 200 on David's.
% fullFieldMicrons = 400;
mirrorFlipY = true; %true for data from Dimitre's rig? def. false for David's 
% mirrorFlipY = true;

autoNormalize = true;
equalize_to_target_histogram =  false;
% numChannels = 1;


if equalize_to_target_histogram
    fs=filesep;
    [ref_file,ref_file_path] = uigetfile({'*.tif';'*.png'},'Please select reference image');
    refimg = imread([ref_file_path,fs,ref_file]);
    refimg_min = double( min(refimg(:)) );
    refimg_max = double(max(refimg(:) ) );
    
    refimg = (double(refimg) - refimg_min)/refimg_max;
    refimg = uint16( refimg * 2^(16-1));
   
    reference_histogram = imhist(refimg);
end

% overlapMethod = 'Max'; %overlaps are calculated as the maximum
overlapMethod = 'Average'; %overlaps are the average of whatever is there
% and the added volume. (not proportional to # regions overlapping the same
% area, weights last added). 

tiffData = dir('*.tif');

% channelacqinfo = cell2mat(regexp(header,'state.acq.imagingChannel\d=[0,1]','match')');
%Let's assume both channels are always on for now. 

if mirrorFlipY
    mF = -1;
else
    mF = 1;
end
%first, determine coordinates and dimensions of combined tile
for m=1:numel(tiffData)
    info = imfinfo(tiffData(m).name);
    header = info(1).ImageDescription;
    
    nchann(m) = length( regexp(header,'state.acq.acquiringChannel\d=1') );
    
    stackpositionX = regexp(header,'state.motor.absXPosition=(-?\d+\.?\d*)','tokens');
    stackpositionX  = str2num(stackpositionX{1}{1})*-1; %mirror flip of X and Y due to reversed motor/image relation
    stackpositionY = regexp(header,'state.motor.absYPosition=(-?\d+\.?\d*)','tokens');
    stackpositionY  = str2num(stackpositionY{1}{1})*mF;
    stackpositionZ = regexp(header,'state.motor.absZPosition=(-?\d+\.?\d*)','tokens');
    stackpositionZ  = str2num(stackpositionZ{1}{1});
    
    zoomFactor = regexp(header,'state.acq.zoomFactor=(\d+\.?\d*)','tokens');
    zoomFactor = str2num(zoomFactor{1}{1});

    zStep = regexp(header,'state.acq.zStepSize=(-*\d+\.?\d*)','tokens');
    zStep = str2num(zStep{1}{1});
    
    zStepA(m) = zStep;
    
%     if zStep<=0
%         error('Negative z step detected!!');
%     end
    
    frames2Average = regexp(header,'state.acq.numberOfFrames=(\d+\.?\d*)','tokens');
    frames2Average = str2num(frames2Average{1}{1});

    trueStacks = numel(info)/(nchann(m)*frames2Average) ;
    
    StackCoordinates(m,:) = [stackpositionX,stackpositionY,stackpositionZ,stackpositionZ+trueStacks*zStep];
    
    if m==1
        fov = fullFieldMicrons/zoomFactor;
        microns_per_pixel = fov/info(1).Width;
    end
    
end
zStep = mode(zStepA);
disp(zStep);
disp(zoomFactor);

numChannels = unique(nchann);

if length( unique(zStepA) ) > 1
    error('multiple step sizes detected!')
end
% zStep = 10;

minPosition = [min(StackCoordinates(:,1)), min(StackCoordinates(:,2)), min(min(StackCoordinates(:,3:4)))];
maxPosition = [max(StackCoordinates(:,1)), max(StackCoordinates(:,2)), max(max(StackCoordinates(:,3:4)))];

EnclosingVolumePx = round([[maxPosition(1:2) - minPosition(1:2)]./microns_per_pixel [maxPosition(:,3)-minPosition(:,3)]/zStep]);


if EnclosingVolumePx(3) < 0
%    EnclosingVolumePx(3) = abs ( EnclosingVolumePx(3) );
   warning('Z Direction is negative! recommend flip z order in original data');
end


% EnclosingVolumePx(1:2) = EnclosingVolumePx(1:2) + info(1).Width + 1; %since the difference is from the center + fudgefactor
EnclosingVolumePx(1:2) = EnclosingVolumePx(1:2) + info(1).Width; 


masterVolume = zeros(EnclosingVolumePx(1),EnclosingVolumePx(2),EnclosingVolumePx(3),'uint16');
mvRef = imref3d(EnclosingVolumePx,[minPosition(2)-fov/2,maxPosition(2)+fov/2],[minPosition(1)-fov/2,maxPosition(1)+fov/2],[minPosition(3),maxPosition(3)]);

%now place each slice in masterVolumeChann
masterVolumeChann1 = masterVolume;
masterVolumeChann2 = masterVolume;

for m=1:numel(tiffData)
info = imfinfo(tiffData(m).name);
header = info(1).ImageDescription;


stackpositionX = regexp(header,'state.motor.absXPosition=(-?\d+\.?\d*)','tokens');
stackpositionX  = str2num(stackpositionX{1}{1})*-1;
stackpositionY = regexp(header,'state.motor.absYPosition=(-?\d+\.?\d*)','tokens');
stackpositionY  = str2num(stackpositionY{1}{1})*mF;
stackpositionZ = regexp(header,'state.motor.absZPosition=(-?\d+\.?\d*)','tokens');
stackpositionZ  = str2num(stackpositionZ{1}{1});

[yP,xP,zP] = worldToIntrinsic(mvRef,stackpositionY,stackpositionX,stackpositionZ);

frames2Average = regexp(header,'state.acq.numberOfFrames=(\d+\.?\d*)','tokens');
frames2Average = str2num(frames2Average{1}{1});

rawImageData = zeros(info(1).Width, info(1).Height, numel(info), 'uint16');

parfor jj =1:numel(info)
   rawImageData(:,:,jj) = imread(tiffData(m).name, jj);
%     rawImageData(:,:,jj) = medfilt2( imread(tiffData(m).name, jj) );
%      rawImageData(:,:,jj) = fliplr( imread(tiffData(m).name, jj) );
end

trueStacks = floor( numel(info)/(numChannels*frames2Average) ); %floor added by Dawnis 09/28/2018. Can deal with incomplete acquisitions. 
averageSubStack = zeros(info(1).Width,info(1).Height, trueStacks,'uint16');

for k=1:numChannels
    channSubStack = rawImageData(:,:,k:numChannels:end);
    for kj = 1:trueStacks
        averageSubStackDBL = mean(channSubStack(:,:,(kj-1)*frames2Average+1:kj*frames2Average),3);
        averageSubStack(:,:,kj) = uint16(averageSubStackDBL);
    end

%     xSpan = [round(xP) - info(1).Width/2+1, round(xP) + info(1).Width/2];
%     ySpan =  [round(yP) - info(1).Height/2+1, round(yP) + info(1).Height/2];
%     zSpan = [round(zP), round(zP) + trueStacks - 1];

    xSpan = [ceil(xP) - info(1).Width/2, ceil(xP) + info(1).Width/2 - 1];
    ySpan =  [ceil(yP) - info(1).Height/2, ceil(yP) + info(1).Height/2 - 1];
    zSpan = [ceil(zP), ceil(zP) + trueStacks - 1];

    if k==1
        masterSubVolume = masterVolumeChann1(xSpan(1):xSpan(2),ySpan(1):ySpan(2),zSpan(1):zSpan(2));
    else
        masterSubVolume = masterVolumeChann2(xSpan(1):xSpan(2),ySpan(1):ySpan(2),zSpan(1):zSpan(2));
    end
    
    if autoNormalize
        maxProject = prctile(averageSubStack(:),99.9);
        minProject = double( median(averageSubStack(:)) );
        
%         yMax = prctile(maxProject(1:end),0.99);
        yMax = max(maxProject(:));
        
        averageSubStack = 2^(16-1)*( [double(averageSubStack) - minProject]./double(yMax) );
        %avoid saturation by normalizing to half range
        
    end
    
    if equalize_to_target_histogram
        parfor jkks =1:size(averageSubStack,3)
            averageSubStack(:,:,jkks) = histeq(averageSubStack(:,:,jkks), reference_histogram);
        end
    end
    
    volumeCombine = double(averageSubStack);
    
    if strcmp(overlapMethod,'Average')
        volumeCombine(masterSubVolume>0) = double( masterSubVolume(masterSubVolume>0) )*0.5 + double( averageSubStack(masterSubVolume>0) )*0.5;
    elseif strcmp(overlapMethod,'Max')
        volumeCombine(masterSubVolume>0) = max( [masterSubVolume(masterSubVolume>0),  averageSubStack(masterSubVolume>0)],[],2);
    end

    if k==1
        masterVolumeChann1(xSpan(1):xSpan(2),ySpan(1):ySpan(2),zSpan(1):zSpan(2)) = uint16( volumeCombine);
    else
        masterVolumeChann2(xSpan(1):xSpan(2),ySpan(1):ySpan(2),zSpan(1):zSpan(2)) = uint16( volumeCombine);
    end
    
    if ~isequal(size(masterVolumeChann1),size(masterVolumeChann2))
        error('Size mismatch between channels. Boo!')
    end
end

end

mkdir('AssembledOutput');
fs = filesep;
metaDataTiff = Tiff(tiffData(1).name,'r');
metaTags = metaDataTiff.getTagNames();
for ch=1:numChannels

    if exist(['AssembledOutput',fs,sprintf('AssembledStackCh%02.0f.tif',ch)],'file')
        delete(['AssembledOutput',fs,sprintf('AssembledStackCh%02.0f.tif',ch)]);
    end
    
    for z=1:size(masterVolume,3)
        tiffObj = Tiff(['AssembledOutput',fs,sprintf('AssembledStackCh%02.0f.tif',ch)],'a');
        
%         tiffObj.setTag('Photometric',metaDataTiff.getTag('Photometric'));
        tiffObj.setTag('Photometric',Tiff.Photometric.LinearRaw);
        tiffObj.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
%         tiffObj.setTag('PlanarConfiguration',metaDataTiff.getTag('PlanarConfiguration'));
        tiffObj.setTag('BitsPerSample',16);
        tiffObj.setTag('SamplesPerPixel',1);
        tiffObj.setTag('Compression',Tiff.Compression.PackBits);
        
        if z==1
            tiffObj.setTag('ImageDescription',metaDataTiff.getTag('ImageDescription'));
        end
%         for fTag=1:numel(metaTags)
%             try
%                 tiffObj.setTag(metaTags{fTag},metaDataTiff.getTag(metaTags{fTag}));
%             catch ME
%                 getReport(ME);
%                 disp(ME);
%             end
% %         tiffObj.setTag('BitsPerSample',BitsPerSample);m
%         end
    
        tiffObj.setTag('ImageWidth',size(masterVolume,2));
        tiffObj.setTag('ImageLength',size(masterVolume,1));
        
        if ch==1
            tiffObj.write(masterVolumeChann1(:,:,z));        
        else
            tiffObj.write(masterVolumeChann2(:,:,z));  
        end
        
        
%         if z==size(masterVolume,3)
            writeDirectory(tiffObj);
            tiffObj.close();
%         end
    end
    
end
    