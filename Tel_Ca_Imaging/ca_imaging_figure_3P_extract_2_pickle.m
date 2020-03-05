clear all, close all
datahome = '/home/dawnis/Dropbox/3PData/Jan30_3P_Tel_CaSeg_TSData/';
[xf,filenames] = xlsread([datahome,'ca_imaging_depths.xlsx']);

image_names_xcl = filenames(2:end,5);
acq_freq = xf(:,3);
img_power = xf(:,1);
depth_microns = xf(:,2);

%filenames for ca_img_data files
datasources = { ...
    '180130_fish1_gcamp6s_telen2_181uW004', '180130_fish1_gcamp6s_OT_161uW009', '180130_fish1_gcamp6s_CB_191uW001';...
    '180130_fish1_gcamp6s_telen2_407uW006', '180130_fish1_gcamp6s_OT_360uW011', '180130_fish1_gcamp6s_CB_191uW2x002';...
    '180130_fish1_gcamp6s_telen2_1324uW008', '180130_fish1_gcamp6s_OT_1600uW012','180130_fish1_gcamp6s_CB_670uW2x003'
    '180130_fish1_gcamp6s_telen2_300uW005',  '180130_fish1_gcamp6s_OT_300uW010',  '';
    '180130_fish1_gcamp6s_telen2_837uW007',  '180130_fish1_gcamp6s_OT_1560uW013', ''};

deep_1Hz = {'180130_fish1_gcamp6s_CB_670uW2x003','180130_fish1_gcamp6s_CB_2130uW1.06Hz2x004'};

directory_parents = {'Tel','OT','CB'};

pickle_path = '/home/dawnis/Data/Data Science/CaTransientCNN/data_3P';

cd(pickle_path);

for region=1:3
   for f = 1:5
       
    if ~isempty(datasources{f,region})
    imgmatch = strcmp(datasources{f,region}, image_names_xcl);
    
    load([datahome,directory_parents{region},'/',datasources{f,region}]);

    mat2np(ca_time_series, [datasources{f,region},'.pkl'], 'float64');
    mat2np(datalabel, [datasources{f,region},'_Label.pkl'], 'uint8');
    
    end
    
   end
end