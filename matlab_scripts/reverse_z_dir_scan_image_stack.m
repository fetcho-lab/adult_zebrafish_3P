%reverse z_direction of a stack.
clear all, close all
[file,path] = uigetfile('*.tif', 'Select stack to reverse the z direction');

cd(path)
info = imfinfo(file);
ScanImage_Header = info(1).ImageDescription;

zStep = regexp(ScanImage_Header,'state.acq.zStepSize=(-*\d+\.?\d*)','tokens');
zxStep = str2num(zStep{1}{1});

replace_start = regexp(ScanImage_Header,'state.acq.zStepSize=','end');
replace_end = regexp(ScanImage_Header,'state.acq.zStepSize=(-*\d+\.?\d*)','end');

if zxStep < 0
    ScanImage_Header(replace_start+1) = [];
else
    ScanImage_Header = [ScanImage_Header(1:replace_start), '-', ScanImage_Header(replace_start+1:end)];
end

line_start = regexp(ScanImage_Header,'state.acq.zStepSize=(-*\d+\.?\d*)','start');
new_end = regexp(ScanImage_Header,'state.acq.zStepSize=(-*\d+\.?\d*)','end');

new_z_step_line = ScanImage_Header(line_start:new_end);
disp(new_z_step_line);
pause;

n_slices = numel(info);

xCnt = 0;

stack = readImage(file);

for x = n_slices:-1:1
    
    xCnt = xCnt + 1;
    
    if xCnt == 1 && exist([file(1:end-4),'reversed_z_step.tif'], 'file')
        delete([file(1:end-4),'reversed_z_step.tif']);
    end
    
    tiffObj = Tiff([file(1:end-4),'reversed_z_step.tif'],'a');

    %         tiffObj.setTag('Photometric',metaDataTiff.getTag('Photometric'));
    tiffObj.setTag('Photometric',Tiff.Photometric.LinearRaw);
    tiffObj.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    %         tiffObj.setTag('PlanarConfiguration',metaDataTiff.getTag('PlanarConfiguration'));
    tiffObj.setTag('BitsPerSample',16);
    tiffObj.setTag('SamplesPerPixel',1);
    tiffObj.setTag('Compression',Tiff.Compression.PackBits);
    tiffObj.setTag('ImageDescription', ScanImage_Header);
    
    tiffObj.setTag('ImageWidth',size(stack,2));
    tiffObj.setTag('ImageLength',size(stack,1));

    slice_to_read = 2*ceil(x/2) - mod(xCnt, 2);
    tiffObj.write(stack(:,:,slice_to_read) );
    
    writeDirectory(tiffObj);
    tiffObj.close();
        
end