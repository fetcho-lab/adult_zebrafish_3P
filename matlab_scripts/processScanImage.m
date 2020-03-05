function [averageSubStack, stackinfo] = processScanImage(filename, channel2get, nAcquiredChannels)
% processed_image = averageSubStackDBL(filename)
% returns an averaged and channel segregated  file along with some
% information from a ScanImage .tif. Channel2Get is the desired channel to
% return in the processed image. 

info = imfinfo(filename);
header = info(1).ImageDescription;

%     stackpositionZ = regexp(header,'state.motor.absZPosition=(-?\d+\.?\d*)','tokens');
stackpositionZ = regexp(header,'state.motor.relZPosition=(-?\d+\.?\d*)','tokens');
stackpositionZ  = str2num(stackpositionZ{1}{1});

zoomFactor = regexp(header,'state.acq.zoomFactor=(\d+\.?\d*)','tokens');
zoomFactor = str2num(zoomFactor{1}{1});

zStep = regexp(header,'state.acq.zStepSize=(-*\d+\.?\d*)','tokens');
zStep = str2num(zStep{1}{1});

stackinfo.zoom = zoomFactor;
stackinfo.z_step = zStep;
stackInfo.z_coordinate = stackpositionZ;

frames2Average = regexp(header,'state.acq.numberOfFrames=(\d+\.?\d*)','tokens');
frames2Average = str2num(frames2Average{1}{1});

parfor jj =1:numel(info)
    rawImageData(:,:,jj) = imread(filename, jj);
    %      rawImageData(:,:,jj) = fliplr( imread(tiffData(m).name, jj) );
end

trueSlices = numel(info)/(nAcquiredChannels*frames2Average);
trueSlices =  floor(trueSlices); %for incomplete stacks

channSubStack = double( rawImageData(:,:,channel2get:nAcquiredChannels:end) );

%channSubStack = channSubStack - min([ median(channSubStack(1,:,1)), median(channSubStack(:,1,1))]); %channel background subtration
%channSubStack(channSubStack < 0) = 0;

%     channSubStack(1:325,:,:) = NaN; %crop edges for some stacks...

averageSubStack = zeros(info(1).Width,info(1).Height,trueSlices);

for kj = 1:trueSlices
    averageSubStackDBL = mean(channSubStack(:,:,(kj-1)*frames2Average+1:kj*frames2Average),3);
%         averageSubStack(:,:,kj) = uint16(averageSubStackDBL);
    averageSubStack(:,:,kj) = averageSubStackDBL;
end

