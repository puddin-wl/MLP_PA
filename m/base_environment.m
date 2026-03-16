% 程序 1：generate_base_scattering.m
% 目标：生成基础介质的散射相差，并记录环境参数 (包含高兼容性 tif 读取)
clear; clc; close all;

%% 1. 物理参数设定
Nx = 1080; Ny = 1080;
lambda = 0.532; % 波长 (um)
p = 1:15;       % Zernike 阶数

% 建立匹配神经网络的归一化对角线坐标系
limit = 1 / sqrt(2);
x = linspace(-limit, limit, Nx);
y = linspace(-limit, limit, Ny);
[X, Y] = meshgrid(x, y);
[theta_full, r_full] = cart2pol(X, Y);

%% 2. 生成介质的散射相位 (Ground Truth 畸变)
% 设定一个较大的畸变系数代表严重的介质散射
base_coeffs_um = gen_zern_coeffs(p, 1.0, 'random1'); 
base_coeffs_rad = (2 * pi / lambda) * base_coeffs_um;
base_waveFront_vec = create_wavefront(p, base_coeffs_rad, r_full(:), theta_full(:));
waveFront_2D_base = reshape(base_waveFront_vec, Ny, Nx);

%% 3. 准备目标图像并补零 (修改为鲁棒读取版)
img_size = 256;

% ！！！请在这里换成你自己的 tif 图片路径！！！
img = imread('volume11.tif'); 

% 3.1 处理多通道问题：如果是 RGB 或多通道图，强制转换为单通道灰度图
if size(img, 3) == 3
    img = rgb2gray(img);
elseif size(img, 3) > 3
    img = img(:, :, 1); 
end

% 3.2 调整图像尺寸到 256x256
img = imresize(img, [img_size, img_size]); 

% 3.3 将图像转为双精度，并强制进行极值归一化 (Min-Max Normalization)
% 这保证了无论是 8位、16位 还是微弱信号，最终像素都会被完美映射到 0~1 之间
img = double(img);
img_min = min(img(:));
img_max = max(img(:));

if img_max > img_min % 防止纯色图片导致除以 0
    img = (img - img_min) / (img_max - img_min);
else
    img = zeros(size(img));
end

% 3.4 创建 1080x1080 的黑画布进行补零
padded_img = zeros(Ny, Nx);
start_idx = (Nx - img_size) / 2 + 1; 
end_idx = start_idx + img_size - 1;  
padded_img(start_idx:end_idx, start_idx:end_idx) = img;

%% 4. 生成一张仅受介质影响的基础模糊图 (用于对比参考)
pupilFun = ones(Ny, Nx) .* exp(1i * waveFront_2D_base);
prf = fftshift(ifft2(ifftshift(pupilFun)));
PSF_base = abs(prf).^2;
PSF_base_norm = PSF_base / sum(PSF_base(:));

O_f = fft2(ifftshift(padded_img));
H_f = fft2(ifftshift(PSF_base_norm));
blurred_padded_base = real(fftshift(ifft2(O_f .* H_f)));
base_blurred_img = blurred_padded_base(start_idx:end_idx, start_idx:end_idx);

%% 5. 保存基础环境供程序 2 使用
save('base_environment.mat', 'base_coeffs_um', 'waveFront_2D_base', ...
     'padded_img', 'Nx', 'Ny', 'start_idx', 'end_idx', ...
     'p', 'lambda', 'r_full', 'theta_full');

% 保存参考图供预览 (为了视觉效果，将浮点矩阵乘回 255 转成 uint8 保存)
imwrite(uint8(base_blurred_img * 255), 'base_scattered_reference.png');

disp('基础散射环境已生成并保存为 base_environment.mat');