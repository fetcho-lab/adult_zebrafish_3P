clear;

selectDirectory = uigetdir(pwd,'Select analysis directory');
cd(selectDirectory);

tiffData = dir('*.tif');
fullField = 400;

for tf = 1:numel(tiffData)

    tiffInfo = imfinfo(tiffData(tf).name);
    
    fprintf('Analyzing %s ... \n ',tiffData(tf).name);
    

    header = tiffInfo(tf).ImageDescription;

    zoomFactor = regexp(header,'state.acq.zoomFactor=(\d+\.?\d*)','tokens');
    zoomFactor = str2num(zoomFactor{1}{1});

    fov = fullField/zoomFactor;

    microns_per_pixel =  fov/tiffInfo(1).Width;
    imgStack = [];
    parfor k=1:numel(tiffInfo)

       imgStack(:,:,k) =  imread(tiffData(tf).name,k) ; 
    end

    imgStack = imgStack(3:end-2,3:end-2,:);
    figure(1),clf
    imagesc(max(imgStack,[],3));

    [x,y] = ginput;

    measureLine = robustfit(x,y);

    lineEndSegments = [ [min(x); max(x)], measureLine(1)+measureLine(2)*[min(x); max(x)]];

    hold on
    plot(lineEndSegments(:,1),lineEndSegments(:,2),'r-');

    figure(2), clf

    clear C
    for k=1:numel(tiffInfo)
    C(k,:) = improfile(imgStack(:,:,k),lineEndSegments(:,1),lineEndSegments(:,2));

    % imagesc(imgStack(:,:,k));
    % pause, clf;

    end

    subplot(1,2,1)
    imagesc(C);

    subplot(1,2,2)
    threshold = median(C(1:end));

    CBinary = medfilt2( C < threshold );
    CBinary = imerode(CBinary,ones(3,3));

    CBinary2 = zeros(size(CBinary),'logical');
    bloodBars = regionprops(CBinary,'Area','PixelList','PixelIdxList');

    timeAxis = [1:size(CBinary,1)]/8.49;
    spaceAxis = [1:size(CBinary,2)]*microns_per_pixel;

    barFit = [];
    for k=1:numel(bloodBars)
       if bloodBars(k).Area > 25
          CBinary2(bloodBars(k).PixelIdxList) = true;

          barFit(:,end+1) = [robustfit(spaceAxis( bloodBars(k).PixelList(:,1) ), timeAxis( bloodBars(k).PixelList(:,2) ) ) ; ...       
                            min(bloodBars(k).PixelList(:,1)); max(bloodBars(k).PixelList(:,1))];
       end
    end


    imagesc(medfilt2( CBinary),'XData',spaceAxis,'YData',timeAxis)
    hold on
    for j=1:size(barFit,2)
       lineConstraints = spaceAxis( [barFit(3,j); barFit(4,j)]);
       plot(lineConstraints , lineConstraints*barFit(2,j) + barFit(1,j),'r-');
    end

    flowRate = median(barFit(2,:))^-1;
    
    fprintf('Measured flow rate was %2.0f microns per second\n',flowRate);
    pause;
end