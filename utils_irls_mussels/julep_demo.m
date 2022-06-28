% Joint reconstruction for multii-band multi-shot diffusion MRI (on Matlab R2018b)
% (c) Erpeng Dai, Stanford University
% Reference: Dai E, et al. Multi-band Multi-shot Diffusion MRI
% Reconstruction with Joint Usage of Structured Low-rank Constraints and Explicit Phase Mapping (JULEP)

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

%% SMS-SENSE and MUSE parameters
lsqr_iter = 200;
lsqr_tol = 1e-3;
lsqr_lambda1 = 5e-5;
lsqr_lambda2 = 1e-2;

par2.lsqr_iter = lsqr_iter;
par2.lsqr_tol = lsqr_tol;
par2.lambda = lsqr_lambda2;

%% load data
filepath='./';
filepath1=[filepath, 'data/'];
savepath=[filepath, 'data/'];
filename=[filepath1, 'MB3_ksp.mat'];
sensname=[filepath1, 'MB3_sensmap.mat'];
saveimg=[savepath, 'img_MB3_julep.mat'];

addpath(genpath('./utils_neatr/'));
addpath(genpath('./utils_irls_mussels/'));
addpath(genpath('./utils_spirit/'));
mkdir(savepath);

load(filename)
par.scale=max(abs(b1ksp2(:)));
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
    %% Reorganize the data, MB slices are concatenated along the RO dimension
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
    
    %% PI for SMS unfolding in only ACS region (VCC is not applicable in this step)
    acs_y=min([64, 2*par.ky_s-par.ky_r, size(b1ksp2, 2)]);% NAV size (default): 64
    acs_x=acs_y*MB;
    sens_acs0=imresize3(sens0, [acs_x, acs_y, num_chan]);
    sens_acs0=sens_acs0./(sos(sens_acs0, 3)+eps);
    %% SMS-SENSE
    if num_sh>1
        ksp1_slc_acs=crop(ksp1_slc, [acs_x, acs_y, num_chan, num_sh]);
    else
        ksp1_slc_acs=crop(ksp1_slc, [acs_x, acs_y, num_chan]);
    end
    par1.sens=sens_acs0;
    par1.lsqr_iter = lsqr_iter;
    par1.lsqr_tol = lsqr_tol;
    par1.N = [acs_x, acs_y];
    par1.num_chan = num_chan;
    par1.lambda = lsqr_lambda1;
    img_sense_acs = single(zeros([acs_x, acs_y, num_sh]));
    for sh = 1:num_sh
        par1.m2d = squeeze(ksp1_slc_acs(:, :, :, sh)~=0);
        
        k = ksp1_slc_acs(:, :, :, sh);
        res = lsqr(@apply_sense_tikc, cat(1, k(:), single(zeros([acs_x*acs_y, 1]))), lsqr_tol, ...
            lsqr_iter, [], [], [], par1);
        img_sense_acs(:, :, sh) = reshape(res, [acs_x, acs_y]);
    end
    
    hanwin0=window2(acs_x, acs_y, @hann);
    hanwin0=single(hanwin0);
    x = ifft2c(fft2c(img_sense_acs).*hanwin0);
    for cnt=1:2
        x = ifft2c(fft2c(x).*hanwin0);
    end
    phasemap=x./abs(x);
    m=abs(mean(img_sense_acs .* conj(phasemap), 3));
    phasemap = phasemap.* (sens_acs0(:,:,1)~=0);
    phasemap=squeeze(phasemap);
    img_sense_acs=m.*phasemap;
    k_recon0=fft2c(img_sense_acs);
    
    %% IRLS MUSSELS
    mask=single(squeeze(sum(abs(ksp1_slc_acs),3))>0);
    [k_recon] =mussels_cs_new(ksp1_slc_acs, sens_acs0, mask, par.ksize, par.niter0, 3, par.lambda, ...
        0, par.tol, k_recon0, 1, vcc_flag);
    
    if vcc_flag
        k_recon=single(k_recon(:, :, 1:par.nSHOT));
    else
        k_recon=single(k_recon);
    end
    img_fista_acs=ifft2c(k_recon);
    
    %% Low-pass filtering
    s2=size(ksp1_slc);
    ksp_fista=fft2c(img_fista_acs);
    ksp_fista=zpad(ksp_fista, [s2(1:2) num_sh]);
    img_fista_full=ifft2c(ksp_fista);
    ksp2_slc=ksp1_slc;
    
    hanwin1=window2(acs_x, acs_y, @hann);
    hanwin1=single(zpad(hanwin1, [s2(1:2)]));
    N=[s2(1)/MB, s2(2)];
    
    %% MUSE Recon
    x = ifft2c(fft2c(img_fista_full).*hanwin1);
    for cnt=1:2
        x = ifft2c(fft2c(x).*hanwin1);
    end
    phasemap=x./abs(x);
    m=abs(mean(img_fista_full .* conj(phasemap), 3));
    phasemap = phasemap.* (sens0(:,:,1)~=0);
    m = m.* (sens0(:,:,1)~=0);
    
    %% VCC
    kspace_joint = single(zeross([N .* [MB,1], num_chan * num_sh]));
    sens = single(zeross([N .* [MB,1], num_chan * num_sh]));
    
    for sh = 1:num_sh
        kspace_coils = squeeze( ksp2_slc(:, :, :, sh) );
        kspace_joint(:,:, 1 + (sh-1) * num_chan : sh * num_chan ) = kspace_coils;
        sens(:,:, 1 + (sh-1) * num_chan : sh * num_chan ) = ...
            dot_mult( sens0, phasemap(:, :, sh));
    end
    
    if vcc_flag
        temp = fft2c(conj(ifft2c( kspace_joint~=0 ))) > 1e-6;
        Kspace_joint = fft2c(conj(ifft2c(kspace_joint))) .* temp;
        Kspace_joint = cat(3, kspace_joint, Kspace_joint);
        par2.sens = cat(3, sens, conj(sens));
    else
        Kspace_joint=kspace_joint;
        par2.sens=sens;
    end
    
    par2.N = N .* [MB,1];
    par2.num_chan = size(par2.sens, 3);
    par2.m2d = (Kspace_joint~=0);
    par2.lambda = lsqr_lambda2;
    
    res = lsqr(@apply_sense_tikc, cat(1, Kspace_joint(:), zeross([prod(N)*MB, 1])), ...
        lsqr_tol, lsqr_iter, [], [], m(:), par2);
    m1 = abs( reshape(res, N.*[MB,1]) );
    
    phasemap=squeeze(phasemap);
    img_muse=m1.*phasemap;
    acs_x=par.kx_r;
    acs_y=par.ky_r;
    sens_acs=repmat(sens0, [1 1 1 num_sh]);
    if num_sh>1
        img_muse=imresize3(img_muse, [acs_x, acs_y, num_sh]);
        ksp1_slc_acs=crop(ksp1_slc, [acs_x, acs_y, num_chan, num_sh]);
    else
        img_muse=imresize(img_muse, [acs_x, acs_y]);
        ksp1_slc_acs=crop(ksp1_slc, [acs_x, acs_y, num_chan]);
    end
    N_acs=size(ksp1_slc_acs);
    
    k_recon0=fft2c(img_muse);
    
    %% IRLS_MUSSELS
    mask=squeeze(sum(abs(ksp1_slc_acs), 3))>0;
    mask=single(mask);
    [k_recon] =mussels_cs_new(ksp1_slc_acs, sens0, mask, par.ksize, par.niter0, par.niter1, par.lambda, ...
        0, par.tol, k_recon0, 0, vcc_flag);
    
    if vcc_flag
        k_recon=single(k_recon(:, :, 1:par.nSHOT));
    else
        k_recon=single(k_recon);
    end
    x_k=ifft2c(k_recon);
    m1_old=m1;
    
    N1=size(k_recon);
    hanwin2=window2(round(N1(1)/2), round(N1(2)/2), @hann);
    hanwin2=single(zpad(hanwin2, [N1(1:2)]));
    phasemap0=ifft2c(k_recon.*hanwin2);
    phasemap0=phasemap0./(abs(phasemap0)+eps);
    m1=abs(mean(x_k.*conj(phasemap0), 3));
    
    sub_img=reshape(m1, [par.kx_r/MB MB par.ky_r]);
    sub_img=permute(sub_img, [1 3 2]);
    fovshift=floor(s2(2)/MB);
    
    ksp1=fftc(sub_img, 2);
    %% CAIPI FOV shift recovery
    phi=2*pi/MB;
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

par.recon_time=toc;
fprintf('Total recon time: %f s\n', par.recon_time);
par.kx_r=par.kx_r/MB;

save(saveimg, 'b1imfinal2', 'par', 'par1', 'par2', '-v7.3')