function PSF = coeffs2PSF(p,coeffs, Sx, pixelSize, lambda, NA, Sz, zStepSize, RI, normFlag)
% calculate PSF (2D or 3D) from Zernike coefficients (OSA/ANSI convention)
%
% Output
%   PSF: 2D or 3D PSF, maximum value normalized to 1
% Input
%   p: a vector of single indexes(OSA/ANSI convention) for Zernike components,
%       elements should be integers(>=0);
%   coeffs: a vector of Zernike coefficients corresponding to p; (unit: um)
%   Sx: image size x and y of the PSF
%   pixelSize: pixel size of the intensity image
%   lambda: wavelength (should have same unit with pixelSize; (unit: um)
%   NA: numerical aperture of objective
%   Sz: image size z of the PSF
%   zStepSize: z pixel size if for 3D PSF 
%   RI: refractive index if for 3D PSF
%   normFlag: 0: no normalization; 1: normalize maximum to 1; [default: 1]

% By: Min Guo
% Apr. 9, 2020
% Oct. 22, 2020: add normFlag choice
Sy = Sx;
if(nargin==6)
    flag3D = 0; %   flag3D: calculate 2D PSF (0) or 3D PSF (1);
    normFlag = 1;
elseif(nargin==7)
    flag3D = 0;
    normFlag = Sz; % turn Sz to nornFlag
elseif(nargin==9)
    flag3D = 1;
    normFlag = 1;
elseif(nargin==10)
    flag3D = 1;
else
    error('coeffs2PSF: incorrect number of input arguments');
end 

% Zernike coefficients: convert lengh unit(um) to phase unit(pi)
length2phase = 2*pi/lambda;
coeffs = length2phase * coeffs;

% ===== 替换开始：生成匹配神经网络的正方形坐标系 =====
% 为了让对角线的最大半径 r 刚好等于 1，x 和 y 的边界应当是 1/sqrt(2)
% 这样在四个角点：r = sqrt((1/sqrt(2))^2 + (1/sqrt(2))^2) = 1
limit = 1 / sqrt(2); 

% 生成等间隔的正方形网格
x = linspace(-limit, limit, Sx);
y = linspace(-limit, limit, Sy);
[X, Y] = meshgrid(x, y);

% 转换为极坐标
[theta_full, r_full] = cart2pol(X, Y);

% 因为没有裁剪，所有的像素都是有效像素
% 我们需要将矩阵展平为列向量传入 create_wavefront
phi_vec = create_wavefront(p, coeffs, r_full(:), theta_full(:)); 

% 将输出的一维向量重新塑形为 2D 正方形矩阵
phi = reshape(phi_vec, Sx, Sy);

% 因为是正方形全通相位，所以 pupilMask 是一个全为 1 的矩阵
pupilMask = ones(Sx, Sy, 'single');
% ===== 替换结束 =====
pupilFun = pupilMask.*exp(1i*phi);
if(flag3D == 0) % 2D PSF at focal plane
    prf = fftshift(ifft2(ifftshift(pupilFun)));
    PSF = abs(prf).^2;
else % 3D PSF
    freSampling = 1/pixelSize; % length^-1
    freSamplingPhase = Sx/freSampling;
    pixelSizePhase = 1/freSamplingPhase;
    % calculate defocus function: dConst
    dConst = zeros(Sx,Sy, 'single');
    Sox = (Sx+1)/2; % Sox == Soy
    for i = 1:Sx
        for j = 1:Sy
            if(pupilMask(i,j)==1)
                rSQ = (i-Sox)^2 + (j-Sox)^2;
                rSQ = rSQ * pixelSizePhase^2;
                dConst(i,j) = sqrt(1-(lambda/RI)^2*rSQ);
            end
        end
    end
    dConst = 2*pi*RI/lambda*dConst;
    % calculate defocus pupil
    PSF = zeros(Sx,Sy,Sz, 'single');
    Soz = (Sz+1)/2;
    for i = 1:Sz
        zPos = (i - Soz) * zStepSize;
        pupilFun2D = pupilFun .*exp(1i*zPos*dConst);
        prf = fftshift(ifft2(ifftshift(pupilFun2D)));
        PSF(:,:,i) = abs(prf).^2;
    end
    % calculate 3D PSF
    
end
% normalization
if(normFlag==1)
    PSF = PSF/max(PSF(:));
end
