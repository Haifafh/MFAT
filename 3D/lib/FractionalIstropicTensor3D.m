function out = FractionalIstropicTensor3D(I,sigmas,tau,tau2,D,spacing,whiteondark)
% calculates the vesselness Fractional anisotropy tensor of a 3D
% input image
%   - Fractional anisotropy tensor equation
%    Reference :
%     1- Hansen, Charles D., and Chris R. Johnson. Visualization handbook. Academic Press, 2011.â€? APA        

% out = FractionalIstropicTensor3D(I, sigmas,spacing,tau ,tau2,D )
% 
% inputs,
%   I : 3D image
%   sigmas : vector of scales on which the vesselness is computed
%   spacing : input image spacing resolution - during hessian matrix 
%       computation, the gaussian filter kernel size in each dimension can 
%       be adjusted to account for different image spacing for different
%       dimensions            
%   tau,tau 2: cutoff thresholding related to eignvlaues.
%   D  : the step size of soultion evolution.
% outputs,
%   out: final vesselness response over scales sigmas
%
% example:
%   out = FractionalIstropicTensor3D(I,[.5:.5:3],1,0.03,0.3,0.27);
%
% Function written by Haifa F. Alhasson , Durham University (Dec 2017)
% Based on code by T. Jerman, University of Ljubljana (October 2014)
%%
spacing = [spacing spacing]; 
verbose = 1;
I(~isfinite(I)) = 0;
%% preprocessing 
I = single(I);
%% Enhancement 
for j = 1:length(sigmas)  
    if verbose
        disp(['Current filter scale (sigma): ' num2str(sigmas(j)) ]);
    end
    %% (1) Eigen-values
    [~, Lambda2, Lambda3] = volumeEigenvalues(I,sigmas(j),spacing,whiteondark);
    % filter response at current scale from RVR
    Lambda3M = Lambda3;
    Lambda3M(Lambda3M<0 & Lambda3M >= tau .* min(Lambda3M(:)))=  tau.* min(Lambda3M(:));
    % New filter response
    Lambda4 = Lambda3;
    Lambda4(Lambda3<0 & Lambda3 >= tau2 .* min(Lambda3(:)))=  tau2.* min(Lambda3(:));
    %% (2) Fractional Anisotropy Tensor equation: 
    % Mean Eigen-value (LambdaMD):
    LambdaMD = (abs(Lambda2)+ abs(Lambda3M) +abs(Lambda4))./3;
    % response at current scale 
    response = sqrt((((abs(Lambda2))-abs(LambdaMD)).^2+(abs((Lambda3M))-abs(LambdaMD)).^2+(abs(Lambda4)-abs(LambdaMD)).^2) ./(((abs(Lambda2))).^2+((abs(Lambda3M))).^2+(abs(Lambda4)).^2));    
    response = sqrt(3./2).*(response);
    response  = imcomplement(response);
    %% (3) restristions:
    x = Lambda3M - Lambda2;
    response(x == max(x(:,:,:))) = 1; 
    response(Lambda3M > x) = 0;
    response(~isfinite(response)) = 0;
    response(Lambda2>=0) = 0;
    response(~isfinite(response)) = 0;
    response(Lambda3M>=0) = 0;   
    response(~isfinite(response)) = 0;   
    %% (4) Update vesselness
    if(j==1)
        vesselness = response;
    else  
        vesselness = vesselness + D .* tanh( abs(response) - D);
        vesselness = max(vesselness,response);
    end
    % Normalize vessleness
    vesselness = min(max(vesselness, 0), 1);
    %% (5) Use vesselness to enhance I for next iteration
%     I = I - tau .* tanh( abs(I-vesselness));
%     I = min(max(I, 0), 1);
    clear Lambda2 Lambda3 Lambda4 LambdaMD 
end
out = vesselness  ./ max(vesselness (:)); % should not be really needed 
out = normalize(out);

function [Lambda1, Lambda2, Lambda3] = volumeEigenvalues(V,sigma,spacing,whiteondark)
% calculates the three eigenvalues for each voxel in a volume

% Calculate 3D hessian
[Hxx, Hyy, Hzz, Hxy, Hxz, Hyz] = Hessian3D(V,sigma,spacing);

% Correct for scaling
c=sigma.^2;
Hxx = c*Hxx; Hxy = c*Hxy;
Hxz = c*Hxz; Hyy = c*Hyy;
Hyz = c*Hyz; Hzz = c*Hzz;

if whiteondark == false
    c=-1;
    Hxx = c*Hxx; Hxy = c*Hxy;
    Hxz = c*Hxz; Hyy = c*Hyy;
    Hyz = c*Hyz; Hzz = c*Hzz;    
end

% reduce computation by computing vesselness only where needed
% S.-F. Yang and C.-H. Cheng, ï¿½Fast computation of Hessian-based
% enhancement filters for medical images,ï¿½ Comput. Meth. Prog. Bio., vol.
% 116, no. 3, pp. 215ï¿½225, 2014.
B1 = - (Hxx + Hyy + Hzz);
B2 = Hxx .* Hyy + Hxx .* Hzz + Hyy .* Hzz - Hxy .* Hxy - Hxz .* Hxz - Hyz .* Hyz;
B3 = Hxx .* Hyz .* Hyz + Hxy .* Hxy .* Hzz + Hxz .* Hyy .* Hxz - Hxx .* Hyy .* Hzz - Hxy .* Hyz .* Hxz - Hxz .* Hxy .* Hyz;

T = ones(size(B1));
T(B1<=0) = 0;
T(B2<=0 & B3 == 0) = 0;
T(B1>0 & B2>0 & B1 .* B2 < B3) = 0;

clear B1 B2 B3;

indeces = find(T==1);

Hxx = Hxx(indeces);
Hyy = Hyy(indeces);
Hzz = Hzz(indeces);
Hxz = Hxz(indeces);
Hyz = Hyz(indeces);
Hxy = Hxy(indeces);

% Calculate eigen values
[Lambda1i,Lambda2i,Lambda3i]=eig3volume(Hxx,Hxy,Hxz,Hyy,Hyz,Hzz);

% Free memory
clear Hxx Hyy Hzz Hxy Hxz Hyz;

Lambda1 = zeros(size(T));
Lambda2 = zeros(size(T));
Lambda3 = zeros(size(T));

Lambda1(indeces) = Lambda1i;
Lambda2(indeces) = Lambda2i;
Lambda3(indeces) = Lambda3i;

% some noise removal
Lambda1(~isfinite(Lambda1)) = 0;
Lambda2(~isfinite(Lambda2)) = 0;
Lambda3(~isfinite(Lambda3)) = 0;

Lambda1(abs(Lambda1) < 1e-4) = 0;
Lambda2(abs(Lambda2) < 1e-4) = 0;
Lambda3(abs(Lambda3) < 1e-4) = 0;


function [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(Volume,Sigma,spacing)
%  This function Hessian3D filters the image with an Gaussian kernel
%  followed by calculation of 2nd order gradients, which aprroximates the
%  2nd order derivatives of the image.
% 
% [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(Volume,Sigma,spacing)
% 
% inputs,
%   I : The image volume, class preferable double or single
%   Sigma : The sigma of the gaussian kernel used. If sigma is zero
%           no gaussian filtering.
%   spacing : input image spacing
%
% outputs,
%   Dxx, Dyy, Dzz, Dxy, Dxz, Dyz: The 2nd derivatives
%
% Function is written by D.Kroon University of Twente (June 2009)
% defaults
if nargin < 2, Sigma = 1; end

if(Sigma>0)
    %F=imbigaussian(Volume,Sigma,0.5);
    F=imgaussian(Volume,Sigma,spacing);
else
    F=Volume;
end

% Create first and second order diferentiations
Dz=gradient3(F,'z');
Dzz=(gradient3(Dz,'z'));
clear Dz;

Dy=gradient3(F,'y');
Dyy=(gradient3(Dy,'y'));
Dyz=(gradient3(Dy,'z'));
clear Dy;

Dx=gradient3(F,'x');
Dxx=(gradient3(Dx,'x'));
Dxy=(gradient3(Dx,'y'));
Dxz=(gradient3(Dx,'z'));
clear Dx;

function D = gradient3(F,option)
% This function does the same as the default matlab "gradient" function
% but with one direction at the time, less cpu and less memory usage.
%
% Example:
%
% Fx = gradient3(F,'x');

[k,l,m] = size(F);
D  = zeros(size(F),class(F)); 

switch lower(option)
case 'x'
    % Take forward differences on left and right edges
    D(1,:,:) = (F(2,:,:) - F(1,:,:));
    D(k,:,:) = (F(k,:,:) - F(k-1,:,:));
    % Take centered differences on interior points
    D(2:k-1,:,:) = (F(3:k,:,:)-F(1:k-2,:,:))/2;
case 'y'
    D(:,1,:) = (F(:,2,:) - F(:,1,:));
    D(:,l,:) = (F(:,l,:) - F(:,l-1,:));
    D(:,2:l-1,:) = (F(:,3:l,:)-F(:,1:l-2,:))/2;
case 'z'
    D(:,:,1) = (F(:,:,2) - F(:,:,1));
    D(:,:,m) = (F(:,:,m) - F(:,:,m-1));
    D(:,:,2:m-1) = (F(:,:,3:m)-F(:,:,1:m-2))/2;
otherwise
    disp('Unknown option')
end
              
function I=imgaussian(I,sigma,spacing,siz)
% IMGAUSSIAN filters an 1D, 2D color/greyscale or 3D image with an 
% Gaussian filter. This function uses for filtering IMFILTER or if 
% compiled the fast  mex code imgaussian.c . Instead of using a 
% multidimensional gaussian kernel, it uses the fact that a Gaussian 
% filter can be separated in 1D gaussian kernels.
%
% J=IMGAUSSIAN(I,SIGMA,SIZE)
%
% inputs,
%   I: 2D input image
%   SIGMA: The sigma used for the Gaussian kernel
%   SPACING: input image spacing
%   SIZ: Kernel size (single value) (default: sigma*6)
% 
% outputs,
%   I: The gaussian filtered image
%

if(~exist('siz','var')), siz=sigma*6; end

if(sigma>0)

    % Filter each dimension with the 1D Gaussian kernels\
    x=-ceil(siz/spacing(1)/2):ceil(siz/spacing(1)/2);
    H = exp(-(x.^2/(2*(sigma/spacing(1))^2)));
    H = H/sum(H(:));    
    Hx=reshape(H,[length(H) 1]);
    
    x=-ceil(siz/spacing(2)/2):ceil(siz/spacing(2)/2);
    H = exp(-(x.^2/(2*(sigma/spacing(2))^2)));
    H = H/sum(H(:));    
    Hy=reshape(H,[1 length(H)]);
    
    I=imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
end

function [Lambda1,Lambda2]=eigvalOfHessian2D(Dxx,Dxy,Dyy)
% This function calculates the eigen values from the
% hessian matrix, sorted by abs value

% Compute the eigenvectors of J, v1 and v2
tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);

% Compute the eigenvalues
mu1 = 0.5*(Dxx + Dyy + tmp);
mu2 = 0.5*(Dxx + Dyy - tmp);

% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
check=abs(mu1)>abs(mu2);

Lambda1=mu1; Lambda1(check)=mu2(check);
Lambda2=mu2; Lambda2(check)=mu1(check);