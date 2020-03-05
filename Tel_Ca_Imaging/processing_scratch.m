close all

[file_select,path] = uigetfile('*.mat','Select Saved Ca Series');
cd(path)
load(file_select);

divide_frequency = true;
smooth_before = true;
frequency = 4.25;

if divide_frequency
    frequency = frequency/2;
    caDS = zeros(size(ca_time_series,1), size(ca_time_series,2)/2);
    for k=2:2:size(ca_time_series,2)
        caDS( :, k/2) = mean( ca_time_series(:,k-1:k),2);
    end
else
    caDS = ca_time_series;
end

if smooth_before
    parfor c = 1:size(caDS,1)
       caDS(c,:) = smooth(caDS(c,:)); 
    end
end

dFF =  pmv_dFF2(caDS, frequency, 60, 50);
maxFF = max(dFF,[],2);
signal_variance = std(dFF,[],2);
[B,I] = sort(maxFF,'descend');

xtick = [1:size(caDS,2)]/frequency;

ccolor = lines(5);

figure(1)
subplot(1,3,1)
img_time_project = mean(imgdata,3);
img_time_project = 2^(16-1)*[double(img_time_project) - min(img_time_project(:))]/max(img_time_project(:));
img_time_project = uint16(img_time_project);
imagesc(img_time_project)
hold on
for k=1:5
   imgblack = zeros(size(imgdata,1), size(imgdata,2),'logical');
   imgblack(spot_segmentation(I(k)).inds) = 1;
   [B,L] = bwboundaries(imgblack);
   plot(B{1}(:,2),B{1}(:,1),'color',ccolor(k,:),'linewidth',1);
end
title('Mean Time Projection');

subplot(1,3,2)
hold on
for k=1:5
plot(xtick, caDS(I(k),:),'color',ccolor(k,:));
end
xlabel('Time (s)');
ylabel('Int');
axis tight
title('Top 5 Fl. Int.');

plotTraces(dFF(I(1:5),:),1,[-8,1],'plotSubWindow',[1,3,3])
subplot(1,3,3)
axis tight
xlabel('Time s')
ylabel('\DeltaF/F');
title('Top 5 \DeltaF/F');