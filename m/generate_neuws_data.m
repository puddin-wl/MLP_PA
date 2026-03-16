clear; clc; close all;

%% 1. 物理参数与系数设定
Nx = 1080; 
Ny = 1080;
lambda = 0.532; % 波长 (um)
p = 1:15;       % 使用前 15 阶 Zernike 模式
aValue = 0.1;   % 系数幅值 (um)

% 生成随机像差系数 (单位: um)
coeffs_um = gen_zern_coeffs(p, aValue, 'random1'); 

%% 2. 坐标系建立 (匹配神经网络对角线逻辑)
x = linspace(-1, 1, Nx);
y = linspace(-1, 1, Ny); 
[X, Y] = meshgrid(x, y);
[theta, r] = cart2pol(X, Y);

% 对角线最大半径归一化为 1
max_r = max(r(:));
r_norm = r / max_r; 

%% 3. 生成相位 (Phase)
% 关键：将微米(um)转换为弧度(rad)，确保物理尺度正确
length2phase = 2 * pi / lambda;
coeffs_rad = length2phase * coeffs_um;

% 生成展平的波前，并重塑为 2D 矩形矩阵
waveFront_vec = create_wavefront(p, coeffs_rad, r_norm(:), theta(:));
waveFront_2D = reshape(waveFront_vec, Ny, Nx);

%% 4. 从相位直接生成 PSF (保证 100% 对应)
% 因为你使用的是正方形全光束，所以振幅掩膜全为 1
pupilFun = 1.0 .* exp(1i * waveFront_2D); 

% 通过逆傅里叶变换计算焦平面的点扩散函数 (PSF)
prf = fftshift(ifft2(ifftshift(pupilFun)));
PSF = abs(prf).^2;

% 归一化 PSF
PSF = PSF / max(PSF(:)); 

%% 5. 可视化验证与保存
figure('Name', 'Matched Phase and PSF', 'Position', [100, 100, 1000, 450]);

% 显示送入神经网络的相位
subplot(1,2,1);
imagesc(waveFront_2D);
colormap(gca, 'jet'); 
colorbar;
title('Saved Phase Pattern (rad)');
axis image; axis off;

% 显示用于卷积图片的 PSF
subplot(1,2,2);
imagesc(PSF);
colormap(gca, 'hot'); 
colorbar;
title('Corresponding 2D PSF');
axis image; axis off;

% 在这里你可以添加代码将 waveFront_2D 和 PSF 导出为 .mat 或图片供 Python 读取
% save('sim_data.mat', 'waveFront_2D', 'PSF');