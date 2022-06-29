% Custom NEATR code based on Dr. Bilgic's open-source code (https://bit.ly/2QgBg9U) to reproduce the
% results in the JULEP paper.
% (1) BM3D is used for denoising; (2) the phase cycling is replaced with a low-pass Hanning window 
% (which demonstrates a better performance on this dataset);
% (c) Erpeng Dai, Stanford University

%% load data
close all
clear
clc

MB=3;
Rpe=1;
num_ech=1;

%% SMS-MUSSELS: Fista iterations
fista = 1;
num_iter = 200;
tol = 0.3;
winSize = [1,1] * 7;
lambda = 2;
%% phase-cycling param
nouteriter = 1;
h = 10;
ninneriter = 50;

%% load data
filepath='./';
filepath1=[filepath, 'data/'];
savepath=[filepath, 'data/'];
filename=[filepath1, 'MB3_ksp.mat'];
sensname=[filepath1, 'MB3_sensmap.mat'];
saveimg=[savepath, 'img_MB3_neatr.mat'];

addpath(genpath('./utils_neatr/'));
mkdir(savepath);

load(filename) % b0ksp2: f0*p0*c0*nex0*s0; b1ksp2: f0*p0*c0*dir0*nex1*s0
par.scale=max(abs(b0ksp2(:)));
b1ksp2=b1ksp2/par.scale;
load(sensname)

if isfield(par, 'blipsign')
    blipsign=par.blipsign;
else
    blipsign=0;
end

[N(1), ~, num_chan, ~, ~, ~]=size(b1ksp2);
N(2)=par.ky_r;

acs_y=min([64, 2*par.ky_s-par.ky_r, size(b1ksp2, 2)]);% NAV size (default): 64
acs_x=acs_y*MB;

%% jvc-sense param
par.lsqr_iter = 200;
par.lsqr_tol = 1e-3; % original value 1e-6 %
par.lambda_L2 = 1e-2; % Tikhonov penalty for Sense

par.kx_r=par.kx_r*MB;
tic

for sl=1:par.nSL
    for nb=1:par.nB
        kspc_colaps=single(zeros([par.kx_r, par.ky_r, par.nCH, par.nSHOT]));
        %% Concatenate MB slices along the RO dimension
        for sh=1:par.nSHOT
            if blipsign>0
                kspc_colaps(1:MB:end, par.ky_r-par.ky_s+(sh-1)*Rpe+1:par.nSHOT*Rpe:end, :, sh)=...
                    b1ksp2(:, (sh-1)*Rpe+1:par.nSHOT*Rpe:end, :, nb, 1, sl);
            else
                kspc_colaps(1:MB:end, (sh-1)*Rpe+1:par.nSHOT*Rpe:par.ky_s, :, sh)=...
                    b1ksp2(:, (sh-1)*Rpe+1:par.nSHOT*Rpe:end, :, nb, 1, sl);
            end
        end
        sens0=sensmap(:, :, :, sl);

        kspc_colaps=permute(kspc_colaps, [1 2 3 5 4]);
        
        kspc_colaps0=kspc_colaps;
        kspc_colaps0(:, par.ky_s+1:par.ky_r, :, :, :)=[];
        kspc_colaps0(:, 1:par.ky_r-par.ky_s, :, :, :)=[];
        N0=N;
        N0(2)=size(kspc_colaps0, 2);
        %% SMS-SENSE
        lsqr_iter = 200;
        lsqr_tol = 1e-3;
        lambda_L1 = 5e-5;
        
        param.sens = imresize3(sens0, [N0 .* [MB 1], num_chan]);
        param.sens = param.sens./sos(param.sens, 3);
        param.N = N0.* [MB,1];
        param.num_chan = num_chan;
        param.lambda = lambda_L1;
        
        img_sense = zeross([N0 .* [MB,1], num_ech, par.nSHOT]);
         
        for sh = 1:par.nSHOT
            param.m2d = sq(kspc_colaps0(:,:,:,1,sh)~=0);
            
            for te = 1:num_ech
                k = kspc_colaps0(:,:,:,te,sh);
                
                res = lsqr(@apply_sense_tikc, cat(1, k(:), zeross([prod(N0) * MB,1])), lsqr_tol, lsqr_iter, [], [], [], param);
                
                img_sense(:,:,te,sh) = reshape(res, N0.* [MB,1]);
            end
        end
        
        %% SMS-MUSSELS: Fista iterations
        Sens = repmat(param.sens, [1,1,1,par.nSHOT]);
        
        cSens = conj(Sens);
        Sens2 = sum(abs(Sens).^2, 3);
        
        img_fista = zeross([MB*N(1), N(2), num_ech, par.nSHOT]);
        
        for ne = 1:num_ech
            k_slice = sq(kspc_colaps0(:,:,:,ne,:));
            
            mask = k_slice~=0;
            
            t_k = 1;
            
            y_k = img_sense(:,:,ne,:);
            
            x_kneg1 = y_k;
            
            for t = 1:num_iter
                % k-space constraint
                im_coils = repmat(y_k, [1,1,num_chan,1]) .* Sens;
                
                im_coils = im_coils + ifft2c3( k_slice - mask .* fft2c3(im_coils));
                
                % sense constraint
                im_pocs = sum(im_coils .* cSens, 3) ./ (eps + Sens2);
                
                % mussel constraint
                A = Im2row( fft2c3(squeeze(im_pocs(:,:,1,:))), winSize );
                [U, S, V] = svd(A, 'econ');
                
                keep = 1:floor(lambda*prod(winSize));
                A = U(:,keep) * S(keep,keep) * V(:,keep)';
                
                k_pocs = Row2im(A, [MB*N0(1), N0(2), par.nSHOT], winSize);
                x_k = ifft2c3(permute(k_pocs, [1,2,4,3]));
                
                if fista
                    t_kplus1 = (1 + sqrt(1 + 4 * t_k^2)) / 2;
                    coef_kneg1 = -(t_k - 1) / t_kplus1;
                    coef_k = (t_kplus1 + t_k - 1) / t_kplus1;
                else
                    coef_k = 1;
                    coef_kneg1 = 0;
                    t_kplus1 = 1;
                end
                
                update = rmse(x_k, x_kneg1);
                disp(['iteration: ', num2str(t), '  update: ', num2str(update), ' %   coef_k: ', num2str(coef_k), '   coef_k-1: ', num2str(coef_kneg1)])
                
                y_kplus1 = coef_k * x_k + coef_kneg1 * x_kneg1;
                t_k = t_kplus1;
                
                y_k = y_kplus1;
                x_kneg1 = y_k;
                
                if update < tol
                    break
                end
            end
            
            img_fista(:,:,ne,:) = imresize(x_k, N.*[MB 1]);
        end
        
        img_sms_mussels=img_fista;
        img_bm3d=[];
        msk_receive=single(abs(sensmap(:, :, 1))>0);
        
        %% bm3d filtering for Sms-Mussels output -> complex valued
        %% uses pre-computed Sms-Mussels recon "img_sms_mussels"
        
        img_bm3d = 0 * img_sms_mussels;
        
        sigma_bm3d = 8;     % filter size
        
        for slc = 1:s(img_sms_mussels,5)
            for sh = 1:s(img_sms_mussels,4)
                disp([slc, sh])
                
                for te = 1:s(img_sms_mussels,3)
                    y = img_sms_mussels(:,:,te,sh,slc);
                    
                    y_re = real(y);
                    min_re = min(y_re(:));
                    y_re = y_re - min_re;
                    
                    y_im = imag(y);
                    min_im = min(y_im(:));
                    y_im = y_im - min_im;
                    
                    scl_re = max(y_re(:));
                    scl_im = max(y_im(:));
                    
                    
                    [~, y_Re] = BM3D(1, y_re / scl_re, sigma_bm3d);
                    [~, y_Im] = BM3D(1, y_im / scl_im, sigma_bm3d);
                    
                    y_est = y_Re * scl_re + min_re + 1i * ( y_Im * scl_im + min_im );
                    
                    img_bm3d(:,:,te,sh,slc) = y_est;
                end
            end
        end
        
        
        img_bm3d = img_bm3d .* repmat(permute(msk_receive,[1,2,3,5,4]),[1,1,1,par.nSHOT,1]);
        
        img_tmp=mean(sos(ifft2c(kspc_colaps0), 3), 5);
        scl_factor=max(img_tmp(:));
        kspc_scale = kspc_colaps / scl_factor;
        
        %% phase cycling & jvc-sense recon -- skipped
        
%         phase-cycling param
%         nouteriter = 1;
%         h = 10;
%         ninneriter = 50;
        
%         C = Identity;
%         M = Identity;
%         P = Identity;
%         S = ESPIRiT( sens0 );
%         
%         
%         ncycles = 16;
%         W = {};
%         
%         lambda_phase = 1e-3;        % L1 wavelet on phase
%         Pp = wave_thresh('db4', 3, lambda_phase);

        m2d = sq(kspc_scale(:, :, :, :)~=0);
        img_jvcsense = zeross([N .* [MB,1], num_ech]);
        
        for te = 1:num_ech
            if sum(abs(img_bm3d(:, :, 1)))==0
                m=mean(abs(img_fista), 4)/scl_factor;
                x=img_fista/scl_factor;
            else
                m=mean(abs(img_bm3d), 4)/scl_factor;
                x=img_bm3d/scl_factor;
            end
            
            
            Phs = sq( angle( x ) ) .* repmat( sens0(:,:,1)~=0, [1,1,par.nSHOT] );
            
            for n = 1:nouteriter
%                 Mm = M * m;
%                 alphap = 1.0 / lipschitz(P) / (max(abs(Mm(:)))^2 + eps) * h;
                
%                 for t = 1:par.nSHOT
%                     disp(['phase cycling for shot: ', num2str(t)])
%                     
%                     y = sq( kspc_scale(:,:,:,te,t) );
%                     
%                     F = p2DFT(m2d(:,:,:,t), [N .* [MB,1], num_chan]);
%                     
%                     p0 = Phs(:,:,t);
%                     p = p0;
%                     
%                     for c = 1:ncycles
%                         p_c = angle( exp(1i * p0) .* exp(1i * (c-1) * 2 * pi /  ncycles) );
%                         
%                         if c == 1
%                             p0_c = p_c;
%                         end
%                         W{c} = p_c - p0_c;
%                     end
%                     
%                     for itp = 1:ninneriter
%                         if isempty(W)
%                             w = 0;
%                         else
%                             ti = randi(length(W));
%                             w = W{ti};
%                         end
%                         
%                         expPp = exp(1j * (P * p));
%                         r = C' * (S' * (F' * (y - (F * (S * (C * (Mm .* expPp)))))));
%                         
%                         if lambda_phase
%                             p = Pp(p + alphap * real(-1j * (P' * (conj(Mm) .* conj(expPp) .* r))) + w, alphap) - w;
%                         else
%                             p = p + alphap * real(-1j * (P' * (conj(Mm) .* conj(expPp) .* r)));
%                         end
%                         
%                     end
%                     
%                     p = p .* (sens0(:,:,1)~=0);
%                     
%                     Phs_result(:,:,t) = p;
%                 end

                Phs_result=Phs;
                img_tmp=exp(1i.*Phs_result);
                hanwin0=window2(acs_x, acs_y, @hann);
                hanwin0=single(zpad(hanwin0, [par.kx_r, par.ky_r]));
                x = ifft2c(fft2c(img_tmp).*hanwin0);
                for cnt=1:2
                    x = ifft2c(fft2c(x).*hanwin0);
                end
                Phs_result=angle(x);
                
                % Enforce real valued recon using conjugate symmetric data
                kspace_joint = zeross([N .* [MB,1], num_chan * par.nSHOT]);
                M2d = zeross([N .* [MB,1], num_chan * par.nSHOT]);
                sens = zeross([N .* [MB,1], num_chan * par.nSHOT]);
                
                for t = 1:par.nSHOT
                    kspace_coils = sq( kspc_scale(:,:,:,te,t) );
                    kspace_joint(:,:, 1 + (t-1) * num_chan : t * num_chan ) = kspace_coils;
                    
                    M2d(:,:, 1 + (t-1) * num_chan : t * num_chan ) = m2d(:,:,:,t);
                    sens(:,:, 1 + (t-1) * num_chan : t * num_chan ) = dot_mult( sens0, exp(1i * Phs_result(:,:,t)) );
                end
                
                % virtual coil data
                temp = fft2c(conj(ifft2c( kspace_joint~=0 ))) > 1e-6;
                Kspace_joint = fft2c(conj(ifft2c(kspace_joint))) .* temp;
                Kspace_joint = cat(3, kspace_joint, Kspace_joint);
                
                param = [];
                param.sens = cat(3, sens, conj(sens));
                
                param.N = N .* [MB,1];
                param.num_chan = size(param.sens,3);
                
                param.m2d = Kspace_joint~=0;
                param.lambda = par.lambda_L2;
                
                res = lsqr(@apply_sense_tikc, cat(1, Kspace_joint(:), zeross([prod(N) * MB, 1])), par.lsqr_tol, par.lsqr_iter, [], [], [], param);
                
                m = abs( reshape(res, N.*[MB,1]) );
            end
            
            img_jvcsense(:,:,te) = m;
        end
        
        sub_img=reshape(img_jvcsense, [par.kx_r/MB, MB, par.ky_r]);
        sub_img=permute(sub_img, [1 3 2]);
        ksp1=fftc(sub_img, 2);

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
        
        sub_img=abs(ifftc(ksp1, 2)) *scl_factor;
        b1imfinal2(:, :, :, nb, sl)=sub_img;
    end
end

par.kx_r=par.kx_r/MB;

save(saveimg, 'b1imfinal2', '-v7.3')