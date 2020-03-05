
%Use this to read stack information from autoAsseembledStack.tif files

[file, path] = uigetfile('*.tif', 'Select stack to read information from')

cd(path)

info = imfinfo(file);
header = info(1).ImageDescription;

stackpositionZ = regexp(header,'state.motor.relZPosition=(-?\d+\.?\d*)','tokens');
stackpositionZ  = str2num(stackpositionZ{1}{1});

zoomFactor = regexp(header,'state.acq.zoomFactor=(\d+\.?\d*)','tokens');
zoomFactor = str2num(zoomFactor{1}{1});

zStep = regexp(header,'state.acq.zStepSize=(-*\d+\.?\d*)','tokens');
zStep = str2num(zStep{1}{1});

 fprintf( 'ZStep: %2.2f, Zoom: %2.0fx, StackZ: %4.2f\n', zStep, zoomFactor, stackpositionZ); 