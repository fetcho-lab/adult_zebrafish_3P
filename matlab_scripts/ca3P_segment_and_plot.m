clear; close all
load F3GCamP6S_3270uW_2xZoom_2.12Hz001.mat

%problem: using background subtracted data!

cell_radius = 5*(200/256)/2;
SE = strel('disk',round(cell_radius),0);

figure(1)
imn = mean(imgdata,3);

groupsize = 10;
imnn=[];
for k=1:size(imgdata,3)/groupsize
   imnn(:,:,k) = mean(imgdata(:,:,(k-1)*groupsize+1:k*groupsize),3); 
end

imagesc(imn);
imnn = max(imnn,[],3);

imnc = localcontrast(uint16(imnn),0.2,0.3);
imnc_ee = imdilate(imerode(imnc,SE),SE);

bwimn = imnc_ee > median(imnc_ee);

a = imnc > 150

A = xcorr2(imerode(imnc,ones(3,3)),double(SE.Neighborhood));

szT = size(SE.Neighborhood);
sz = szT(1);
asz = size(A); asz = asz(1);

Av = A(floor(sz)/2:asz-floor(sz)/2,floor(sz)/2:asz-floor(sz)/2);
Avo = Av;

%place cells until the area of the cells falls under threshold
mask = zeros(size(Av));
cellcount = 0;
while max(Av(:)) > 200
    mxAv = max(Av(:));
    [cx,cy] = find(Av == mxAv);
  
    rangeHy = [cy-floor(sz/2):cy+floor(sz/2)];
    rangeVx = [cx-floor(sz/2):cx+floor(sz/2)];
    
    a = 1:5;
    b = 1:5;
    
    a(rangeHy<1) = [];
    b(rangeVx<1) = [];
    a(rangeHy>size(Av,1) ) = [];
    b(rangeVx>size(Av,1) ) = [];
    
    sen = SE.Neighborhood(b,a);

    rangeHy(rangeHy<1) = [];
    rangeVx(rangeVx<1) = [];
    
    rangeHy(rangeHy>size(Av,1) ) = [];
    rangeVx(rangeVx>size(Av,1) ) = [];
    
    currentMaskValue = mask(rangeVx,rangeHy);
    
    if sum(currentMaskValue(sen))  > 0
        nptt = currentMaskValue;
        nptt(~currentMaskValue & sen) = cellcount+1;
        mask(rangeVx,rangeHy) = nptt;
        Av(rangeVx,rangeHy) = 0;
        
    else
        mask(rangeVx,rangeHy) = double(sen)*(cellcount+1);
        Av(rangeVx,rangeHy) = 0;
    end
    
    cellcount = cellcount+1;
%     imagesc(Av), pause, clf;
end

cellProperties = regionprops(mask,imnn,'Area','PixelIdxList','MeanIntensity','Centroid');
AreaThreshold = 11;
intThreshold = 250;

Area = [cellProperties.Area];
Int = [cellProperties.MeanIntensity];

valid = [Area > AreaThreshold] &  [Int > intThreshold];

cellsValid = cellProperties(valid);

lblImg = zeros(size(Avo));

for m=1:numel(cellsValid)
    lblImg(cellsValid(m).PixelIdxList) = m;
end

labelImg = label2rgb(lblImg);

figure(1),clf
subplot(2,1,1)
imagesc(imnn);
axis equal
centers = vertcat(cellsValid.Centroid);
hold on, viscircles(centers,ones(size(centers,1),1)*2.5,'LineWidth',1,'EnhanceVisibility',false);

subplot(2,1,2)
imshow(labelImg);

%can add manual selection like in the decay analysis. 

save Segmentation3270uW_001 cellsValid imnn Av

cellF =  zeros(numel(cellsValid),800);

parfor k=1:800
    
    slice_o_time = zeros(numel(cellsValid),1);
    imgg = imgdata(:,:,k);
    
    for c=1:numel(cellsValid)
         
        slice_o_time(c) = mean( imgg (cellsValid(c).PixelIdxList ) );
    end
    
    cellF(:,k) = slice_o_time;
    
end

dFF = pmv_dFF2(cellF,2.12,150,50);

for k=1:44
    group = [(k-1)*10 + 1 : k*10];
    dFFrw = dFF(group,:);
    
%     dFFsm = rowfun(@(x) smooth(x), dFFrw,'OutputFormat','uniform');

    parfor j=1:10
        dFFsm(j,:) = smooth(dFFrw(j,:));
    end
    
    figure(1), cla
    plotTraces(dFFsm,1,[-10 2],'plotSubWindow',[2,1,1]);
    
    subplot(2,1,2)
    FRw = cellF(group,:);
    
    parfor j=1:10
      Fsm(j,:) = smooth(FRw(j,:));
    end
    
    plot(bsxfun(@plus, Fsm, -[0:9]'*200)');
    
    pause, clf
end

