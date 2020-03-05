function deltaFF = pmv_dFF2_3P(timeseries,stackfrequency,windowsize,basepercent)
%deltaFF = pmv_dFF2(timeseries,stackfrequency,windowsize,basepercent)
%Computes a moving baseline deltaF/F for timeseries dat. windowsize is the
%time in seconds over which the baseline will be computed and basepercent
%is the percentile of the data in the window which will be used as the
%baseline. 'parallel' specifies that parallel processing should be used.

windowStacks = round(stackfrequency*windowsize);
[nCells,nStacks] = size(timeseries);

if basepercent < 1
    warning('Base percent is 50 for 50%');
end
    
indexGrid = zeros(nStacks,windowStacks);
for g=1:size(indexGrid,1)
   indexGrid(g,:) = g-floor(windowStacks/2)+1:g+ceil(windowStacks/2 );
end

indexGrid(indexGrid<1) = NaN;
indexGrid(indexGrid>nStacks) = NaN;

indexGridPAD = isnan(indexGrid);
rows2iterate = sum(indexGridPAD,2) > 0;

nonpadded_iG = indexGrid(~rows2iterate,:);
toPad = find(rows2iterate);


% nonpadded_iG = indexGrid(ceil(windowStacks/2):nStacks-ceil(windowStacks/2),:);
fbaseline = zeros(1,nStacks);

maxFreq = stackfrequency/5; %attenuate high frequency peaks associated with movement artifacts
nOrder = round(3*stackfrequency); 
NyQuist = stackfrequency/2;
[B,A] = butter(nOrder,maxFreq/NyQuist,'low');
    
%     if exist('toParallel','var') && strcmp(toParallel,'parallel')
        parfor c=1:nCells
            cellF = timeseries(c,:);        
            cellF_lp = filtfilt(B,A,cellF);
            deltaFF(c,:) = doOneCell(indexGrid,nonpadded_iG,toPad,cellF_lp,basepercent,rows2iterate,fbaseline);
        end
%     else
%         for c=1:nCells
%             cellF = timeseries(c,:);        
%             cellF_lp = filtfilt(B,A,cellF);
%             deltaFF(c,:) = doOneCell(indexGrid,nonpadded_iG,toPad,cellF,cellF_lp,basepercent,rows2iterate,fbaseline);
%         end
%     end
    
    
end

function dFFc = doOneCell(indexGrid,nonpadded_iG,toPad,lowPass,basepercent,rows2iterate,fbaseline)
%     cellSliding = cellF(nonpadded_iG);
    cellSliding = lowPass(nonpadded_iG);

    baselineMiddle = prctile(cellSliding',basepercent);
    fbaseline(~rows2iterate) = baselineMiddle;
    
%     for k = 1:length(toPad)
%         gIdx = toPad(k);
%         if gIdx < length(cellF/2)
%             fbaseline(gIdx) = baselineMiddle(1);
%         else
%             fbaseline(gIdx) = baselineMiddle(end);
%         end
%     %     validIdx = ~isnan(indexGrid(gIdx,:));
%     %     fbaseline(gIdx) = prctile(cellF(indexGrid(gIdx,validIdx)),basepercent);
%     end
    windowStacks = size(indexGrid,2);
    nStacks = length(lowPass);
    
    for k=1:floor(windowStacks/2)-1
        validIdx = ~isnan(indexGrid(k,:));
        baselineStart(k) = prctile(lowPass(indexGrid(k,validIdx)),basepercent);
    end
    
    for k=nStacks-ceil(windowStacks/2)+1:nStacks
        validIdx = ~isnan(indexGrid(k,:));
        baselineEnd(k - (nStacks-ceil(windowStacks/2))) = prctile(lowPass(indexGrid(k,validIdx)),basepercent);
    end
    
    fbaseline = zeros(1,nStacks);
    fbaseline(rows2iterate) = [baselineStart baselineEnd];
    fbaseline(~rows2iterate) = baselineMiddle;
    
    dFFc = [lowPass-fbaseline]./fbaseline;
    dFFc = [dFFc - median(dFFc)]; %since a percentile != 0.5 will cause a shift in the baseline
end