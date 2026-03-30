% 程序 2：generate_neuws_dataset.m
% 目标：在底层散射环境中叠加 SLM 随机相位，生成 100 组网络训练数据
clear; clc; close all;

%% 1. 加载底层散射环境 (程序 1 生成的结果)
load('base_environment.mat'); 
% 这会载入: base_coeffs_um, waveFront_2D_base, padded_img, Nx, Ny, start_idx, end_idx, p, lambda, r_full, theta_full

%% 2. 实验参数设置
num_samples = 100;     % 生成 100 组数据
slm_aValue = 0.5;      % SLM 随机调制的幅值 (um)，可以比底层畸变(1.0)略小
output_dir = 'dataset_neuws'; % 保存数据的文件夹名称

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

disp(['开始生成 ', num2str(num_samples), ' 组数据集...']);

%% 3. 循环生成数据
for idx = 1:num_samples
    
    % 3.1 生成第 idx 个 SLM 的随机相位
    slm_coeffs_um = gen_zern_coeffs(p, slm_aValue, 'random1'); 
    slm_coeffs_rad = (2 * pi / lambda) * slm_coeffs_um;
    
    % 注意：强烈建议你在 create_wavefront 函数里强制使用 'NOLL' 排序，以匹配 Python 源码！
    % 如果 create_wavefront 还不支持传参，请确保底层调用的 zernfun 是 Noll 排序。
    slm_waveFront_vec = create_wavefront(p, slm_coeffs_rad, r_full(:), theta_full(:), 'norm', 'NOLL');
    waveFront_2D_slm = reshape(slm_waveFront_vec, Ny, Nx);
    
    % 3.2 物理光学叠加：总相位 = 组织散射相位 + SLM 相位
    % (对应网络前向传播中的 H_aber * H_slm)
    total_waveFront = waveFront_2D_base + waveFront_2D_slm;
    
    % 3.3 计算叠加后的系统总 PSF
    pupilFun_total = ones(Ny, Nx) .* exp(1i * total_waveFront);
    prf_total = fftshift(ifft2(ifftshift(pupilFun_total)));
    PSF_total = abs(prf_total).^2;
    PSF_total_norm = PSF_total / sum(PSF_total(:)); % 能量守恒归一化
    
    % 3.4 卷积成像 (频域加速)
    O_f = fft2(ifftshift(padded_img));
    H_f_total = fft2(ifftshift(PSF_total_norm));
    blurred_padded_total = real(fftshift(ifft2(O_f .* H_f_total)));
    
    % 3.5 截取中心有效区域 (256x256 探测器/重构区域)
    final_captured_img = blurred_padded_total(start_idx:end_idx, start_idx:end_idx);
    
    % 3.6 将 SLM 相位缩小至 256x256 (为了送入网络)
    % 网络期望的输入分辨率是 256x256，而不是原始的 1080x1080
    proj_sim = imresize(waveFront_2D_slm, [256, 256]);
    imsdata = final_captured_img;
    
    % 3.7 保存为 NeuWS 网络要求的 .mat 格式和变量名
    slm_filename = fullfile(output_dir, sprintf('SLM_sim%d.mat', idx));
    img_filename = fullfile(output_dir, sprintf('SLM_raw%d.mat', idx));
    
    save(slm_filename, 'proj_sim');
    save(img_filename, 'imsdata');
    
    % 打印进度
    if mod(idx, 10) == 0
        fprintf('已生成 %d / %d 组数据\n', idx, num_samples);
    end
end

disp('数据集生成完毕！可以送给 Python 端训练了。');

%% 4. 可视化最后一组的结果做个检查
figure('Name', 'Sample Checking', 'Position', [200, 200, 1200, 400]);

subplot(1,3,1);
imagesc(waveFront_2D_base); colormap(gca, 'jet'); colorbar; axis image off;
title('Base Tissue Aberration');

subplot(1,3,2);
imagesc(waveFront_2D_slm); colormap(gca, 'jet'); colorbar; axis image off;
title(['SLM Modulation #', num2str(num_samples)]);

subplot(1,3,3);
imagesc(imsdata); colormap(gca, 'gray'); colorbar; axis image off;
title(['Captured Speckle Image #', num2str(num_samples)]);
