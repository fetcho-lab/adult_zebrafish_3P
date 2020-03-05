clear all, close all
fs=filesep;
%cd ('C:\Users\Joe\Dropbox\Jan30_3P_Tel_CaSeg_TSData') 
cd ('/home/dawnis/Dropbox/3PData/Jan30_3P_Tel_CaSeg_TSData');
[xf,filenames] = xlsread('ca_imaging_depths.xlsx');

image_names_xcl = filenames(2:end,5);
acq_freq = xf(:,3);
img_power = xf(:,1);
depth_microns = xf(:,2);

%filenames for ca_img_data files
datasources = { ...
    '180130_fish1_gcamp6s_telen2_181uW004', '180130_fish1_gcamp6s_OT_161uW009', '180130_fish1_gcamp6s_CB_191uW001';...
    '180130_fish1_gcamp6s_telen2_407uW006', '180130_fish1_gcamp6s_OT_360uW011', '180130_fish1_gcamp6s_CB_191uW2x002';...
    '180130_fish1_gcamp6s_telen2_1324uW008', '180130_fish1_gcamp6s_OT_1600uW012','180130_fish1_gcamp6s_CB_670uW2x003'};

deep_1Hz = {'180130_fish1_gcamp6s_CB_670uW2x003','180130_fish1_gcamp6s_CB_2130uW1.06Hz2x004'};

directory_parents = {'Tel','OT','CB'};


original_data_directory = '/home/dawnis/Data/3PImaging/1300nm/Jan30_2018_GCamP6S/30';

load cell_display_order;

for region=1:3
   for f = 1:3
    imgmatch = strcmp(datasources{f,region}, image_names_xcl);
    
    load([directory_parents{region},fs,datasources{f,region}]);
%   

    npx = arrayfun(@(x) numel(x.inds), spot_segmentation);
    npx_selected = npx( logical(datalabel));
    
    original_data_file = [original_data_directory,filesep,datasources{f,region},'.tif'];
    tifData = imread(original_data_file);
    oData_header = imfinfo(original_data_file);
    
    scan_image_header = oData_header(1).ImageDescription;
    chann3_auto = regexp(scan_image_header,'state.acq.pmtOffsetAutoSubtractChannel3=0');
    
    if isempty(chann3_auto)
        error('Auto bg subtract on detected for file');
    end
    
    instrument_black = mean( mean( tifData(:, 1:2)) );
    
    frame_median = zeros( 1, numel(oData_header));
    parfor k=1:numel(oData_header)
        img = imread(original_data_file, k);
        frame_median(k) = median(img(:));
    end
    
    if strcmp(datasources{f,region}, '180130_fish1_gcamp6s_OT_161uW009')
        frame_median(1065:end) = []; %hand correction for manually truncated data
    end
    correction_factor = double( instrument_black - frame_median );
    correction_factor_half_freq = 0.5*correction_factor(1:2:end) + 0.5*correction_factor(2:2:end) ;
    
    figure(2)
    subplot(3,3,(f-1)*3+region)
    freq = acq_freq(imgmatch);
    good_series = ca_time_series(logical(datalabel),:);
    good_series_half_freq = 0.5* (good_series(:,1:2:end) + good_series(:,2:2:end));
    good_series_half_freq = smoothdata(good_series_half_freq,2,'movmean',5);
%     photon_count = pmv_photon_count2(good_series_half_freq, freq/2, 60, 50);

    good_series_half_freq_corrected = bsxfun(@plus, good_series_half_freq, correction_factor_half_freq);
    photon_count = freq*bsxfun(@times, good_series_half_freq_corrected, npx_selected')/880;
    
    title(directory_parents{region});
    ylabel([ sprintf('Depth %3.0f',depth_microns(imgmatch)), ' \mum']);
    
    photCount_Save = [directory_parents{region},'_',num2str(depth_microns(imgmatch)),'_Depth_photonCount.mat'];
    save(photCount_Save,'photon_count','-v7.3');
    imagesc(photon_count)
%     caxis([-0.05 2]);
%     ylabel(sprintf('Depth %3.0f',depth_microns(imgmatch)));
    title([sprintf('%s: Depth %3.0f',directory_parents{region},depth_microns(imgmatch)), '\mum'] );
    ylabel('Cell#');
    set(gca,'XTick',[0:100:800],'XTicklabel',round( [0:100:800]/(freq/2) ));
    xlabel('Time (s)');
    colormap jet
    colorbar;
    
    figure(3)
    good_idx = find(datalabel);
%     maxSignal = max(photon_count,[],2);
%     [sMax,sIndex] = sort(maxSignal,'descend');
    
    sIndex = sIndex_storage{region}{f};
    
    nsignals = min(size(photon_count,1), 10);
%     fH = plotTraces(photon_count(sIndex(1:nsignals),:),1,[-12 1],'plotSubWindow',[3,3,(f-1)*3+region]);
    subplot(3,3,(f-1)*3 + region);
    fH = plot(photon_count(sIndex(1:nsignals),:)');
    set(gca,'XTick',[0:100:800],'XTicklabel',round( [0:100:800]/(freq/2) ));
    title(directory_parents{region});
    title([sprintf('%s: Depth %3.0f',directory_parents{region},depth_microns(imgmatch)), '\mum'] );
    ylabel('Photon Flux (counts per s)');
%     ylabel([ sprintf('Depth %3.0f',depth_microns(imgmatch)), ' \mum']);
    xlabel('Time (s)');
    axis tight
    
%     
%     figure(5+region)
%     subplot(1,3,f)
    
    figure(5)
    subplot(3,3,(f-1)*3 + region)
    
    img_time_project = mean(imgdata,3);
    img_time_project = 2^(16-1)*[double(img_time_project) - min(img_time_project(:))]/max(img_time_project(:));
    img_time_project = uint16(img_time_project);
    imagesc(img_time_project)
    colormap(jet);
    title(directory_parents{region});
    ylabel([ sprintf('Depth %3.0f',depth_microns(imgmatch)), ' \mum']);
    set(gca,'XTick',[],'YTick',[]);
    axis equal
    axis tight
    hold on
    
%     for k=1:size(photon_count,1)
%        imgblack = zeros(size(imgdata,1), size(imgdata,2),'logical');
%        imgblack(spot_segmentation(good_idx(k)).inds) = 1;
%        [B,L] = bwboundaries(imgblack);
%        plot(B{1}(:,2),B{1}(:,1),'color','r','linewidth',0.2);
% 
% %         cctr = double( spot_segmentation(good_idx(k)).center );
% %         text(cctr(1),cctr(2),num2str(k),'color','r');
%     end
    
    ccolor = lines(nsignals);
    gS_Idx = good_idx(sIndex);
    for k=1:nsignals
           imgblack = zeros(size(imgdata,1), size(imgdata,2),'logical');
           imgblack(spot_segmentation(gS_Idx(k)).inds) = 1;
           [B,L] = bwboundaries(imgblack);
           plot(B{1}(:,2),B{1}(:,1),'color',ccolor(k,:),'linewidth',0.5);
    end
    
    
    %now circle the cells which are plotted in figure 3 in figure 1
    
    
   end
end