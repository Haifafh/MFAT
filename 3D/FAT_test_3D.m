%% Clear
clc; clear; close all; 
%% Patch
addpath('..\3D\images\');
addpath('..\3D\lib\');
addpath('..\3D/results\');
%% Folders
LFAT = 'results\LFAT\';
PFAT = 'results\PFAT\';
ext1 = '.png';
ext2 = '.mat';
foldername = '\images\';
%% Add Filename
filename = '1lsm';
ext = '.tif';
%% Load image
I = load3D([foldername filename ext]); 
%% Normalise original image
I = (I - min(I(:)))/(max(I(:))-min(I(:))); 
%I = imcomplement(I);
%I = rgb2gray(I); 
%% Parameters setting
sigmas = [.5:.5:2];
spacing =.7;
whiteondark= true;
tau = 0.02;
tau2 = 0.2;
D = 0.18;
%% Proposed Method (Eign values based version)
iout1 = FractionalIstropicTensor3D(I,sigmas,tau,tau2,D,spacing,whiteondark);
iout1 = normalize(iout1);
clear Hxx Hyx Hxy Hyy
clear Lambda2 Lambda3 Lambda4 response x
%% Proposed Method (probability based version)
iout2 = ProbabiliticFractionalIstropicTensor3D(I,sigmas,tau,tau2,D,spacing,whiteondark);
iout2 = normalize(iout2);
clear Hxx Hyx Hxy Hyy
clear Lambda2 Lambda3 Lambda4 response x
%% Display output
figure,subplot(1,3,1),imagesc(max(I,[],3)), title(' Original image');colormap gray; axis equal; axis off;
subplot(1,3,2),imagesc(max(iout1,[],3)),title('LFAT'); colormap gray; axis equal; axis off;
subplot(1,3,3),imagesc(max(iout2,[],3)),title('PFAT');colormap gray; axis equal; axis off;
%% Save images
imwrite(max(iout1,[],3),(fullfile(LFAT,[filename ext1])),'XResolution',300);
imwrite(max(iout2,[],3),(fullfile(PFAT,[filename ext1])),'XResolution',300);
imwrite(max(I,[],3),(fullfile(Original,[filename ext1])),'XResolution',300);
