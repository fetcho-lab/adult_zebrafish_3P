%script for power curve analysis with different aperature sizes in red on
%Feb01 and 21, 2018

cd('/home/dawnis/Data/3PImaging/1700nm/Feb21_2018/PowerCurve_OT');

fData = dir('*.tif');

datatable = [];

for k=1:numel(fData)
   file_name = fData(k).name;
   
   recorded_aperature_cell = regexp(file_name,'([0-9]+)mm','tokens');
   recorded_aperature = str2double( recorded_aperature_cell{1}{1} );
   
   recorded_power_cell = regexp(file_name,'([0-9]+\.[0-9])mW','tokens');
   
   if ~isempty(recorded_power_cell)
       recorded_power = str2double(recorded_power_cell{1}{1});
   else
       recorded_power = 0;
   end
   
   region = regexp(file_name,'telen','match');
   if ~isempty(region)
       imageregion = 1;
   else
       imageregion = 0;
   end
   
   stackdata = processScanImage(file_name, 1, 2);
%    stackdata  = stackdata - median(stackdata(:));
%    stackdata(stackdata<0) = 0;
   brightest = prctile(stackdata(:), 99);
   
%    power_corrected = recorded_power * [1-exp(-2 * 21.5 / recorded_aperature) ] ;
    power_corrected = PMeter_ap_calc(recorded_power, recorded_aperature);
%    log_normalized_signal = log( brightest/ (power_corrected^3) );
%    log_normalized_signal_no_correct = log( brightest/ (recorded_power^3) );
   
                        %1                  2                   3               4           5            6
   datatable(k,:) = [recorded_aperature, recorded_power, power_corrected, imageregion, brightest]; % log_normalized_signal, log_normalized_signal_no_correct];
   
%    imagesc(max(stackdata,[],3) );
%    pause; clf;
end


laseroff = find( datatable(:,2) == 0 );
% datatable(laseroff, 6) = NaN; %log(datatable(laseroff, 5) );
% mean_laseroff = mean(datatable(laseroff,5));
laseroff_val = datatable(laseroff(1),5);
signal_l = datatable(:,5) - laseroff_val;
datatable(:,6) = log(signal_l ./ datatable(:,3).^3 );
datatable(:,7) = log(signal_l ./ datatable(:,2).^3 );
datatable(laseroff,:) = [];

aperature_sizes = unique(datatable(:,1) );
aperature_cmap = hsv( numel(aperature_sizes) );

power_cmap = hsv(5);
[histn,edges,bin] = histcounts(datatable(:,3), 5);

figure(1), clf

subplot(1,2,1)

for m=[0,1]
    region = datatable(:,4) == m;
    
    datasubtable = datatable(region,:);
    bincode = bin(region);
    
    hold on
%     subplot(1,2,m+1);

    y = scatter(datasubtable(:,1), datasubtable(:,6),50,power_cmap(bincode,:),'filled');
    if m== 0
        h = plot(datasubtable(:,1), datasubtable(:,6), 'ko', 'linewidth', 2);
    end
%     for a = 1:numel(aperature_sizes)
%         relevant = datasubtable(:,1) == aperature_sizes(a);
%         
%         power = datasubtable(relevant,3);
%         [power, sI] = sort(power);
%         
%         signal = datasubtable(relevant,6);
%         signal = signal(sI);
%         
%         if power(1) == 0
%             x = plot(power(1),signal(1),'rx','linewidth',2);
%             power(1) = [];
%             signal(1) = [];
%         end
%         if m == 0
%             h = plot(power, signal, 'o-', 'color', aperature_cmap(a,:), 'markerfacecolor', aperature_cmap(a,:));
%         else
%             y = plot(power, signal, 'o--', 'color', aperature_cmap(a,:), 'markerfacecolor', aperature_cmap(a,:));
%         end
%     end
% 
% if m==0
%     legend ( [h,x], {'OT','laseroff'});
% else
%     legend ( [y], {'Tel'});
% end
% 
% xlabel('Corrected Power (mW)');
% ylabel('log (int/power^3)');
% 
% if m==1
    colormap(power_cmap)
    cb = colorbar;
    set(cb,'XTick',[0:0.2:1], 'XTickLabel', edges);
    ylabel(cb, 'Power (mW)');
% end

% axis([0 50 -4 8]);
end
% legend ( [h,y,x], {'OT','Tel','laseroff'});
% 
title('Corrected Power');
legend([h,y], {'Optic Tectum', 'Telencephalon'})
xlabel('Diameter of Aperature (mm)');
ylabel('log (int/power^3)');

% colormap(aperature_cmap)
% cb = colorbar;
% set(cb,'XTick',[0.125:0.25:1], 'XTickLabel', aperature_sizes);
% ylabel(cb, 'Aperature Size (mm)');


[histn,edges,bin] = histcounts(datatable(:,2), 5);
subplot(1,2,2), hold on
for m=[0,1]
    region = datatable(:,4) == m;
    
    datasubtable = datatable(region,:);
    bincode = bin(region);
    
    y = scatter(datasubtable(:,1), datasubtable(:,7),50,power_cmap(bincode,:),'filled');
    if m==0
        h = plot(datasubtable(:,1), datasubtable(:,7), 'ko', 'linewidth', 2);
    end

    colormap(power_cmap)
    cb = colorbar;
    set(cb,'XTick',[0:0.2:1], 'XTickLabel', edges);
    ylabel(cb, 'Power (mW)');

end
title('Recorded Power');
xlabel('Diameter of Aperature (mm)');
ylabel('log (int/power^3)');
legend([h,y], {'Optic Tectum', 'Telencephalon'})