%% Clear
clc; clear; close all; 
%% Adding Paths
addpath('../2D/images/');
addpath('../2D/lib/');
addpath('../2D/results/');
%% Folders
LFAT = 'results/LFAT/';
PFAT = 'results/PFAT/';
ext1 = '.png';
ext2 = '.mat';
foldername = './images/';
%% Add Filename
filename = 'sampleImage';
ext = '.jpg';
%% Load image
I = imread([foldername filename ext]); 
%% Normalize the image
I = normalize(I);
I = mat2gray(imcomplement(I(:,:,2)));
%% Parameters setting
sigmas = [1:1:5];
spacing = .7; whiteondark= true;
tau = 0.05; tau2 = 0.25; D = 0.45;
%% Proposed Method (Eign values based version)
iout1 = FractionalIstropicTensor(I,sigmas,tau,tau2,D,spacing,whiteondark);
iout1 = normalize(iout1);
clear Hxx Hyx Hxy Hyy
clear Lambda2 Lambda3 Lambda4 response  x

%% Proposed Method (probability based version)
iout2 = ProbabiliticFractionalIstropicTensor(I,sigmas,tau,tau2,D,spacing,whiteondark);
iout2 = normalize(iout2);
clear Hxx Hyx Hxy Hyy
clear Lambda2 Lambda3 Lambda4 response x

%% Display output
figure,subplot(1,3,1),imagesc(I), title('Original image');
colormap gray; axis equal; axis off;
subplot(1,3,2),imagesc(iout1),title('LFAT'); 
colormap gray; axis equal; axis off;
subplot(1,3,3),imagesc(iout2),title('PFAT'); 
colormap gray; axis equal; axis off;

%% Save and Write
%Save as 'png' image
imwrite(iout1,(fullfile(LFAT,[filename ext1])),'XResolution',300);
imwrite(iout2,(fullfile(PFAT,[filename ext1])),'XResolution',300);