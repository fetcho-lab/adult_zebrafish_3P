clear all, close all
cd ('/home/dawnis/Dropbox/3PData/Jan30_3P_Tel_CaSeg_TSData')
[xf,filenames] = xlsread('ca_imaging_depths.xlsx');

brainRegion = filenames(2:end,1);
dataName = filenames(2:end,5);

for k=9:numel(dataName)
    
    fullpath = [brainRegion{k},'/',dataName{k},'.mat'];
    load(fullpath);
    
    time_vector = [1:size(ca_time_series,2)]/xf(k,3);
    datalabel = zeros(size(ca_time_series,1),1);
    
    for m=1:size(ca_time_series,1)
        
        signalmax = max(ca_time_series(m,:));
        plot(time_vector, smooth(ca_time_series(m,:)));
        ylim([0 max([400,signalmax])]);
      
        datalabel(m) = input('1 for good, 0 for bad: ');
        
        
    end
    
    save(fullpath, 'datalabel', '-append');
    fprintf('Saving labels for %s\n', fullpath);
end