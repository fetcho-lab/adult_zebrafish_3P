%decay length analysis on raw folder
clear; close all;

analfolder = uigetdir(pwd,'Select Analysis Folder');
cd(analfolder);

wv = input('Is imaging at 1300nm or 1700nm? ');
%%%%%%%%%%%%%%%%%%%%%%%%%%user parameters

if wv == 1300
   fullField=400; %downstairs
   powerScaleFactor = 60;
elseif wv == 1700
    fullField=230; %upstairs
    powerScaleFactor = 1;
end

% minArea = 12;
minArea = 16;


numChannels = 2;

nDilate = 2;
nErode = 1;
deStrel = ones(3,3);

stdAboveMean = 1.96;  %threshold of binary above mean
% stdAboveMean = 1.2;  %threshold of binary above mean

zOffSets = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%updated to read powers from filenames...
% [NMBR,TXTS] = xlsread('stack_information.xlsx');
% datafiles = TXTS(2:end,1);
% pwrSeries = NMBR(:,1);
% zMax = NMBR(1,2);
% 
% 
% 
% if zOffSets
%    zOff = NMBR(:,3);
% end

tifFiles = dir('*.tif');

datafiles = {tifFiles.name};

AA = cellfun(@(x) regexp(x,'(\d+)[u,m]W','tokens'), datafiles, 'UniformOutput',false);
pwrSeries = cell2mat( cellfun(@(x) str2num(x{1}{1}), AA, 'UniformOutput',false) );

zMax = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



decayData = []; %each row has, average intensity, depth, and power as a triple. 

for tf = 1:length(datafiles)
    
    %raw image extraction and averaging
    info = imfinfo(datafiles{tf});
    header = info(1).ImageDescription;
    
%     stackpositionZ = regexp(header,'state.motor.absZPosition=(-?\d+\.?\d*)','tokens');
    stackpositionZ = regexp(header,'state.motor.relZPosition=(-?\d+\.?\d*)','tokens');
    stackpositionZ  = str2num(stackpositionZ{1}{1});
    
    if zOffSets
        stackpositionZ + zOff(tf);
    end
    
    zoomFactor = regexp(header,'state.acq.zoomFactor=(\d+\.?\d*)','tokens');
    zoomFactor = str2num(zoomFactor{1}{1});
    
    zStep = regexp(header,'state.acq.zStepSize=(-*\d+\.?\d*)','tokens');
    zStep = str2num(zStep{1}{1});
    
    frames2Average = regexp(header,'state.acq.numberOfFrames=(\d+\.?\d*)','tokens');
    frames2Average = str2num(frames2Average{1}{1});
    
    fov = fullField/zoomFactor;
    
    parfor jj =1:numel(info)
        rawImageData(:,:,jj) = imread(datafiles{tf}, jj);
        %      rawImageData(:,:,jj) = fliplr( imread(tiffData(m).name, jj) );
    end
    
    trueSlices = numel(info)/(numChannels*frames2Average);
    trueSlices =  floor(trueSlices); %for incomplete stacks
    channSubStack = double( rawImageData(:,:,1:2:end) );
    
    channSubStack = channSubStack - min([ median(channSubStack(1,:,1)), median(channSubStack(:,1,1))]); %channel background subtration
    channSubStack(channSubStack < 0) = 0;
    
%     channSubStack(1:325,:,:) = NaN; %crop edges for some stacks...
    
    averageSubStack = zeros(info(1).Width,info(1).Height,trueSlices);
    
	for kj = 1:trueSlices
        averageSubStackDBL = mean(channSubStack(:,:,(kj-1)*frames2Average+1:kj*frames2Average),3);
%         averageSubStack(:,:,kj) = uint16(averageSubStackDBL);
        averageSubStack(:,:,kj) = averageSubStackDBL;
    end
    
    %image processing and value extraction
    sliceDepth = []; slicePower = []; pixelIntensityAverageRw = []; pixelBright = [];
    for kj=1:size(averageSubStack,3)
        imgData = double(averageSubStack(:,:,kj) );
%         medFiltImg = medfilt2(imgData);
%         
%         dilateImg = medFiltImg;
% %         dilateImg = imgData;
%         for dl = 1:nDilate
%             dilateImg = imdilate(dilateImg,deStrel );
%             
%         end
%         
%         erodeImg = dilateImg;
%         for er=1:nErode
%            erodeImg = imerode(erodeImg,deStrel ); 
%         end
% %         dilate_erode_img = imerode( imdilate(medFiltImg, ones(5,5) ), ones(5,5) );
%         
%         dilate_erode_img = erodeImg;
        inputResponse = 'R';
        
        while ~isempty(inputResponse)
            pixelHi = imgData > nanmean(imgData(1:end))  + 3.291*nanstd(imgData(1:end)) ; %should get top 0.1% of pixels in each frame;
            pixelBright(kj) = nanmean(imgData(pixelHi));
            
            if strcmp(inputResponse,'R')
                bvBinary = imgData >  nanmedian(imgData(1:end) ) + 	stdAboveMean*nanstd(imgData(1:end) );
                
                bvBinary = medfilt2(bvBinary);
                
                for dl = 1:nDilate
                    dilateImg = imdilate(bvBinary,deStrel );
                    
                end
                
                erodeImg = dilateImg;
                for er=1:nErode
                   erodeImg = imerode(erodeImg,deStrel ); 
                end
                bvBinary = erodeImg;
                
            end

            statStructure = regionprops(bvBinary,'Area','PixelIdxList','Centroid');

            deleteStruct = [];
            for k=1:numel(statStructure)
                if statStructure(k).Area * (fov/info(1).Width)^2 <  minArea %get rid of areas less than minArea square microns
                    bvBinary(statStructure(k).PixelIdxList) = false;
                    
                    deleteStruct = [deleteStruct, k];
                end
                
            end
            statStructure(deleteStruct) = [];
            
            pixelIntensityAverageRw(kj) = mean(imgData(bvBinary));
            sliceDepth(kj) = zMax - stackpositionZ - zStep*(kj-1);
            slicePower(kj) = pwrSeries(tf);

            figure(1),clf
            subplot(1,2,1)
            imshow(imgData);
            caxis([min(imgData(1:end)),max(imgData(1:end))] );
            title(sprintf('Depth %2.0f um, Power %2.2f mW', sliceDepth(kj), slicePower(kj)/powerScaleFactor) );

            subplot(1,2,2)
            imshow(bvBinary);
            hold on,
            for kjk = 1:numel(statStructure)
                text(statStructure(kjk).Centroid(1),statStructure(kjk).Centroid(2),num2str(kjk),'Color','g','FontSize',10);
            end

            inputResponse = input('press RETURN to accept or r to reject segmentation, R to reset, or # to remove patch: ','s');

            if strcmp(inputResponse,'r')
                pixelIntensityAverageRw(kj) = NaN;
                break;
            elseif strcmp(inputResponse,'#')
               subplot(1,2,2);
               Button=1;
               disp('Left click to remove, right click to end');
               while Button ~=3
                   [x,y,Button] = ginput(1);
                   
                   if Button == 3
                       break;
                   elseif isempty(x)
                       disp('Stop hitting return, dummy');
                       continue; 
                   end
                   
                   centroid_list = vertcat(statStructure.Centroid);
                   distance_to_click = pdist2([x,y],centroid_list);
                   
                   [min_distance, toRemove] = min(distance_to_click);
                   
                   if min_distance < 30
                        
                       bvBinary(statStructure(toRemove).PixelIdxList) = false;
                       statStructure(toRemove) = [];
                       
                        subplot(1,2,2),cla
                        imshow(bvBinary);
                        hold on,
                        for kjk = 1:numel(statStructure)
                            text(statStructure(kjk).Centroid(1),statStructure(kjk).Centroid(2),num2str(kjk),'Color','g','FontSize',10);
                        end
                   end
                   
               end
            end
        end
        
    end
    
    decayData = [decayData; pixelIntensityAverageRw' pixelBright' sliceDepth' slicePower' ones(trueSlices,1)*tf];
    
end
%minimum blood vessel area is 15 um^2
save decayAnalysis decayData

figure,
numstacks = numel( unique(decayData(:,5)) );
cmap = hsv(numstacks);
% p0 = min(decayData(:,4));
hold on
for tj=1:numstacks
stackIdx = decayData(:,5) == tj;
% dataLog = log (decayData(stackIdx,1)./ [(decayData(stackIdx,4)/p0).^3]);
dataLog = log (decayData(stackIdx,1)./ [decayData(stackIdx,4).^3]);
plot(decayData(stackIdx,3),dataLog,'ko','markerfacecolor',cmap(tj,:));

% stackFit(tj,:) = robustfit(decayData(stackIdx,3),log(decayData(stackIdx,1)./ [(decayData(stackIdx,4)/p0).^3]));
    try
        stackFit(tj,:) = robustfit(decayData(stackIdx,3),log(decayData(stackIdx,1)./ [decayData(stackIdx,4).^3]));
    catch
        stackFit(tj,:) = NaN;
    end
end

pxH=plot(decayData(:,3),log(decayData(:,2)./ (decayData(:,4).^3)),'b.');
ylabel('Log Av. Int. / I^3');
xlabel('Depth (\mum)');
title('Decay Curve');

maxDepth = max(decayData(:,3));
minDepth = min(decayData(:,3));

curveFit = robustfit(decayData(:,3), log(decayData(:,1)./[decayData(:,4).^3]));
fH = plot([minDepth; maxDepth], [minDepth,maxDepth]*curveFit(2)+curveFit(1),'k-');
legend([fH,pxH], {sprintf('Log S =  %2.3f + %2.3f*L',curveFit(1), curveFit(2)),'Top 0.1% Pixels'} );


