%%
%Get File and Put it in Imaris
[file_to_get, path] = uigetfile('*.mat', 'Select registered calcium file');
cd (path)

load(file_to_get);
% img_time_project = max(imgdata,[],3);
img_time_project = mean(imgdata,3);
img_time_project = 2^(16-1)*[double(img_time_project) - min(img_time_project(:))]/max(img_time_project(:));
img_time_project = uint16(img_time_project);
imagesc(img_time_project)

im=GetImaris;

j = im.GetFactory.CreateDataSet; 
j.Create(Imaris.tType.eTypeUInt16,size(imgdata,1),...
         size(imgdata,2),1,...
         1, 1); 
j.SetUnit('um');
j.SetExtendMaxX(400); j.SetExtendMaxY(400); j.SetExtendMaxZ(5);
j.SetExtendMinX(0); j.SetExtendMinY(0); j.SetExtendMinZ(0);
im.SetDataSet(j);

PutVolume(img_time_project, j, 0, 0);
%%
%Run Segmentation
segment_HuCH2BGCamP_ims_3P


%%
%Extract spot data and calcium time series
im=GetImaris;
fs=filesep;
im.FileSave([path, fs, file_to_get(1:end-4),'.ims'], 'writer="Imaris5"');
clear spot_segmentation
sp_objs = CheckObjects(im, 'Spots');
obj_to_get = strcmp('caSegmentation',sp_objs(:,1));


segmented_spots = sp_objs{obj_to_get,2};
spots_centers = segmented_spots.GetPositionsXYZ;
spots_radii = segmented_spots.GetRadiiXYZ;

%use Imaris to get pixels
colDim = size(img_time_project,1);
rowDim = size(img_time_project,2);
[dX,dY] = meshgrid(0:colDim-1,0:rowDim-1); %switched 07/22/2016
vdX = dX*Sc(1,1)+Sc(1,1)/2; %convert each pixel on the meshgrid to a voxel center point
vdY = dY*Sc(2,2)+Sc(2,2)/2;

maxEnclosure = max(spots_radii(:,1:3),[],1);
mxEnclPix = ceil(maxEnclosure./diag(Sc)');
mxEnclPix = mxEnclPix(1:2);
sample_box_dims =2*mxEnclPix+2*[2,2]; %include a buffer in pixel-space
    
for k=1:size(spots_centers,1)
    ell = diag(spots_radii(k,1:2).^2)^-1;
    
    aPOffX = round(spots_centers(k,1)/Sc(1,1) - sample_box_dims(1)/2):round(spots_centers(k,1)/Sc(1,1) + sample_box_dims(1)/2); aPOffX(aPOffX<0) = []; aPOffX(aPOffX>colDim-1) = [];
    aPOffY = round(spots_centers(k,2)/Sc(2,2) - sample_box_dims(2)/2):round(spots_centers(k,2)/Sc(2,2) + sample_box_dims(2)/2); aPOffY(aPOffY<0) = []; aPOffY(aPOffY>rowDim-1) = [];
    subIndxX = dX(aPOffY+1,aPOffX+1)+1; subIndxY = dY(aPOffY+1,aPOffX+1)+1;
    
    total_elements_n = numel(aPOffX) * numel(aPOffY);

%     linear_idx = sub2ind([rowDim,colDim], reshape(subIndxY,total_elements_n,1),...
%                                        reshape(subIndxX,total_elements_n,1));
    linear_idx = sub2ind([rowDim,colDim], reshape(subIndxX,total_elements_n,1),...
                                       reshape(subIndxY,total_elements_n,1));

  
    subVol_voxel_centers = [vdX(linear_idx),vdY(linear_idx)];
    x_xc = zeros(size(subVol_voxel_centers));
    
    for m =1:2
        x_xc(:,m) = subVol_voxel_centers(:,m)-spots_centers(k,mod(m,2)+1);
    end
    
    dist2Ellipse = sum([x_xc*ell] .* x_xc,2);
    enclosed_solution = dist2Ellipse < 1;
    
   spot_segmentation(k).inds = linear_idx(enclosed_solution);
   spot_segmentation(k).center = spots_centers(k,1:2);
end

% % % use segmentation script to get pixels
% segment_centers = vertcat(cell_info.center)*Sc(1,1);
% for k=1:size(spots_centers,1)
%     spot_match = pdist2(spots_centers(k,1:2), segment_centers);
%     match =  spot_match == 0;
%     
%    spot_segmentation(k) = cell_info(match);
% end

%%
clear imgmax_mean_int
clear imaris_mean_int
imaris_mean_int = GrabIMStatistics(segmented_spots,'Spots', 'Intensity Mean');
imgmaxcopy = img_time_project;
for k=1:numel(spot_segmentation)
   imgmax_mean_int(k) = mean(img_time_project(spot_segmentation(k).inds)); 
   imgmaxcopy(spot_segmentation(k).inds) = 0;
end
figure(1),clf
subplot(1,2,1)
plot(imaris_mean_int(:,3),imgmax_mean_int,'ko','markerfacecolor','k');
xlim = get(gca, 'xlim');
hold on; plot(xlim,xlim,'r-');
xlabel('Imaris Intensity Mean');
ylabel('Matlab Intensity Mean');
title('Spot position sanity check');


subplot(1,2,2)
imagesc(img_time_project);
hold on
for k=1:numel(spot_segmentation)
   imgblack = zeros(size(img_time_project),'logical');
   imgblack(spot_segmentation(k).inds) = 1;
   [B,L] = bwboundaries(imgblack);
   plot(B{1}(:,2),B{1}(:,1),'w','linewidth',1);
end
title('Segmentation Ported Into Matlab');

% subplot(1,3,3)
% imagesc(imgmaxcopy);

%%
ca_time_series = zeros(numel(spot_segmentation), size(imgdata,3));

parfor t=1:size(imgdata,3)
   time_slice = imgdata(:,:,t);
   ca_time_series_column = zeros(numel(spot_segmentation),1);
   
   for k=1:numel(spot_segmentation)
      ca_time_series_column(k) = mean( time_slice(spot_segmentation(k).inds) );
   end
   
   ca_time_series(:,t) = ca_time_series_column;
end
% figure(1), clf
% plot(ca_time_series');

save(file_to_get,'imgdata','ca_time_series','spot_segmentation','spots_centers');

%%
% figure(2)
fm = bsxfun(@minus, ca_time_series, mean(ca_time_series,2));
fm_half = zeros(size(fm,1), size(fm,2)/2);
even = [1:2:size(fm,2)];
odd = [2:2:size(fm,2)];

for m=1:size(fm,1)
   fm_half(m,:) = 0.5* (fm(m,even) + fm(m,odd) );
   fm_half(m,:) = smooth(fm_half(m,:)) ;
    
end
% subplot(1,2,1),
plotTraces(fm, 100, [-2000 100]);
