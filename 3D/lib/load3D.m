function [I,fname] = load3D(fname,pname)
%% open 3D tiff
if nargin==0
    [fname, pname] = uigetfile({'*.tif;*.tiff','TIFF Image Files (*.tiff, *.tif)'},'Select the TIFF image file');
elseif nargin==1
    pname = sprintf('%s',pwd);
end
input = sprintf('%s%s',pname, fname);
info = imfinfo(input);
num_images = numel(info);
I = zeros(info(1).Height, info(1).Width, num_images);
for i = 1:num_images
    I(:,:,i)=imread(input,i,'Info',info);
end
clear i info num_images

end