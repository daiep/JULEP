% Implementation of IRLS MUSSELS for comparison (on Matlab R2018b)
% (c) Erpeng Dai, Stanford University

close all
clear
clc

%% Input parameters
MB=3;
Rpe=1;
vcc_flag=1;

%% IRLS MUSSELS parameters
ksize=[4 4];
niter0=5;
niter1=40;
lambda=3e-4;
tol=5e-3;

%% load data
filepath='./';
filepath1=[filepath, 'data/'];
savepath=[filepath, 'data/'];
filename=[filepath1, 'MB3_ksp.mat'];
sensname=[filepath1, 'MB3_sensmap.mat'];
saveimg=[savepath, 'img_MB3_julep.mat'];

addpath(genpath('./utils_neatr/'));
addpath(genpath('./utils_irls_mussels/'));
mkdir(savepath);

load(filename) % b0ksp2: f0*p0*c0*nex0*s0; b1ksp2: f0*p0*c0*dir0*nex1*s0
par.scale=max(abs(b0ksp2(:)));
b1ksp2=b1ksp2/par.scale;

load(sensname)
par.kx_r=par.kx_r*MB;
par.vcc_flag=vcc_flag;
par.ksize=ksize;
par.niter0=niter0;
par.niter1=niter1;
par.lambda=lambda;
par.tol=tol;

if isfield(par, 'blipsign')
    blipsign=par.blipsign;
else
    blipsign=0;
end

b1imfinal2=single(zeros([par.kx_r/MB, par.ky_r, MB, par.nB, par.nSL]));
tic;

for sl=1:par.nSL
    for nb=1:par.nB   
    num_chan=par.nCH;
    num_sh=par.nSHOT;
    %% Concatenate MB slices along the RO dimension
    ksp1_slc=single(zeros([par.kx_r, par.ky_r, par.nCH, par.nSHOT]));
    for sh=1:par.nSHOT
        if blipsign>0
            ksp1_slc(1:MB:end, par.ky_r-par.ky_s+(sh-1)*Rpe+1:par.nSHOT*Rpe:end, :, sh)=...
                b1ksp2(:, (sh-1)*Rpe+1:par.nSHOT*Rpe:end, :, nb, 1, sl);
        else
            ksp1_slc(1:MB:end, (sh-1)*Rpe+1:par.nSHOT*Rpe:par.ky_s, :, sh)=...
                b1ksp2(:, (sh-1)*Rpe+1:par.nSHOT*Rpe:end, :, nb, 1, sl);
        end
    end
    sens0=sensmap(:, :, :, sl);
    
    ksp1_slc_acs=ksp1_slc;
    k_recon0=[];
    %% IRLS_MUSSELS
    mask=squeeze(sum(abs(ksp1_slc_acs), 3))>0;
    mask=single(mask);
    [k_recon] =mussels_cs_new(ksp1_slc_acs, sens0, mask, par.ksize, par.niter0, par.niter1, par.lambda, ...
        0, par.tol, k_recon0, 1, vcc_flag);
    
    if vcc_flag
        k_recon=single(k_recon(:, :, 1:par.nSHOT));
    else
        k_recon=single(k_recon);
    end
    x_k=ifft2c(k_recon);
    
    N1=size(k_recon);
    hanwin2=window2(round(N1(1)/2), round(N1(2)/2), @hann);
    hanwin2=single(zpad(hanwin2, [N1(1:2)]));
    phasemap0=ifft2c(k_recon.*hanwin2);
    phasemap0=phasemap0./(abs(phasemap0)+eps);
    m1=abs(mean(x_k.*conj(phasemap0), 3));
    
    sub_img=reshape(m1, [par.kx_r/MB MB par.ky_r]);
    sub_img=permute(sub_img, [1 3 2]);
    fovshift=floor(par.ky_r/MB);
    
    ksp1=fftc(sub_img, 2);
    %% CAIPI FOV shift recovery
    phi=2*pi/MB; % or -2*pi/MB depending on CAIPI implementation
    for mb=1:MB % slice index
        for kz=1:MB % kz index
            kz_idx=mod(floor(MB/2)+(kz-1), MB);
            kz_off=-(kz_idx-(MB-1)/2);
            phi1=kz_off*phi*((mb-1)-floor(MB/2));
            if par.blipsign>0
                ksp1(:, kz:MB:end, mb)=...
                    ksp1(:, kz:MB:end, mb)*exp(-1j*phi1);
            else
                ksp1(:, end-kz+1:-MB:1, mb)=...
                    ksp1(:, end-kz+1:-MB:1, mb)*exp(-1j*phi1);
            end
        end
    end
    
    sub_img=abs(ifftc(ksp1, 2));
    b1imfinal2(:, :, :, nb, sl)=sub_img;
    end
end

par.kx_r=par.kx_r/MB;

save(saveimg, 'b1imfinal2', '-v7.3')