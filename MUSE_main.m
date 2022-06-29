% Implementation of MUSE for comparison (on Matlab R2018b)
% (c) Erpeng Dai, Stanford University

close all
clear
clc

%% Input parameters
MB=3;
Rpe=1;
vcc_flag=1;

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
saveimg=[savepath, 'img_MB3_muse.mat'];

addpath(genpath('./utils_neatr/'));
addpath(genpath('./utils_irls_mussels/'));
mkdir(savepath);

load(filename) % b0ksp2: f0*p0*c0*nex0*s0; b1ksp2: f0*p0*c0*dir0*nex1*s0
par.scale=max(abs(b0ksp2(:)));
b1ksp2=b1ksp2/par.scale;

load(sensname)
par.kx_r=par.kx_r*MB;
par.vcc_flag=vcc_flag;

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
    
    %% Low-pass filtering
    s2=size(ksp1_slc);
    ksp_irls=fft2c(img_sense_acs);
    ksp_irls=zpad(ksp_irls, [s2(1:2) num_sh]);
    img_irls_full=ifft2c(ksp_irls);
    ksp2_slc=ksp1_slc;
    
    hanwin1=window2(acs_x, acs_y, @hann);
    hanwin1=single(zpad(hanwin1, [s2(1:2)]));
    N=[s2(1)/MB, s2(2)];
    
    %% MUSE Recon
    x = ifft2c(fft2c(img_irls_full).*hanwin1);
    for cnt=1:2
        x = ifft2c(fft2c(x).*hanwin1);
    end
    phasemap=x./abs(x);
    m=abs(mean(img_irls_full .* conj(phasemap), 3));
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
    
    sub_img=reshape(m1, [par.kx_r/MB MB par.ky_r]);
    sub_img=permute(sub_img, [1 3 2]);
    fovshift=floor(s2(2)/MB);
    
    ksp1=fftc(sub_img, 2);
    %% CAIPI FOV shift recovery
    phi=2*pi/MB;% or -2*pi/MB depending on CAIPI implementation
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