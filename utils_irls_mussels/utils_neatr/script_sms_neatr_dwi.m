%--------------------------------------------------------------------------
%% add path 
%--------------------------------------------------------------------------

addpath(genpath('/cluster/kawin/Congyu/CODE/msEPI/NEATR_Siemens/'))   

% addpath /cluster/kawin/berkin/Matlab_Code_New/NEATR/

a_setpath_v2;
 


file_path = '/autofs/space/lucifer_001/users/data/2018_11_15_Bay4_InVivo_diffusion_4subjects_NEATR/subject2/';
file_name = 'meas_MID00303_FID66811_msEPI_1mm_PAT9_9shot.dat';

f_DiffData = [file_path, file_name];         
% shot_accl = [1:9]        % indices of shots to use in Mussels


for diff_dir = 0:6
    shot_accl = [1:9] + 9 * diff_dir       % indices of shots to use in Mussels
    

    %load EPI_data
    out = mapVBVD(f_DiffData);

    n=numel(out);
    x=out{n};
    y= x.image{:,:,:,:,shot_accl,:};    % data

    clear meas
    
    meas.data(:,:,:,1,:,1,1,1,:,1,:)=y;
    clear x y

    meas.data = permute(meas.data,[1 3 2 10 8 7 9 11 4 5 6]);

    meas.data_phascor1d = out{n}.phasecor(); %navigator for data
    meas.data_phascor1d = permute(meas.data_phascor1d,[1 3 2 10 8 7 9 11 4 5 6]);

    s = size(meas.data_phascor1d);
    s(2) = 3; s(11) = 1;
    test = zeros(s);
    test(:,1,:,:,:,:,:,1,:,:) = meas.data_phascor1d(:,:,:,:,:,:,:,1,:,:,1);
    test(:,2,:,:,:,:,:,2,:,:) = meas.data_phascor1d(:,:,:,:,:,:,:,2,:,:,1);
    test(:,3,:,:,:,:,:,1,:,:) = meas.data_phascor1d(:,:,:,:,:,:,:,1,:,:,2);

    meas.data_phascor1d = test;
    meas.data = mrir_image_slice_deinterleave(meas.data);  %deinterleave
    meas.data_phascor1d = mrir_image_slice_deinterleave(meas.data_phascor1d);

    %load patrefscan
    meas.patrefscan = out{n}.refscan(); %in-plane reference scan
    meas.patrefscan = permute(meas.patrefscan,[1 3 2 10 8 7 9 11 4 5 6]);
    meas.patrefscan_phascor = out{n}.refscanPC(); % navigator for patrefscan
    meas.patrefscan_phascor = permute(meas.patrefscan_phascor,[1 3 2 10 8 7 9 11 4 5 6]);

    s = size(meas.patrefscan_phascor);
    s(2) = 3; s(11) = 1;
    test = zeros(s);
    test(:,1,:,:,:,:,:,1,:,:) = meas.patrefscan_phascor(:,:,:,:,:,:,:,1,:,:,1);
    test(:,2,:,:,:,:,:,2,:,:) = meas.patrefscan_phascor(:,:,:,:,:,:,:,2,:,:,1);
    test(:,3,:,:,:,:,:,1,:,:) = meas.patrefscan_phascor(:,:,:,:,:,:,:,1,:,:,2);
    meas.patrefscan_phascor = test;

    %deinterleave
    meas.patrefscan = mrir_image_slice_deinterleave(meas.patrefscan);
    meas.patrefscan_phascor = mrir_image_slice_deinterleave(meas.patrefscan_phascor);

    [meas.prot meas.evp] = read_meas_prot(f_DiffData);
    clear out s test



    %--------------------------------------------------------------------------
    % ghost correct
    %--------------------------------------------------------------------------

    PAT_acqMethod='std'; %'std','FLEET'


    AccZ = 1;        %SMS factor
    AccY = 9;         %PAT factor

    offset_ky = 0;
    del_ky = 0:AccY-1 
    num_shot=length(shot_accl);
    startLine=2;

    nSlice = size(meas.data,10)*AccZ;
    nCoil = size(meas.data,3);
    nLine = size(meas.data,1);
    nColumn = size(meas.data,2);
    nColumnData = round(size(meas.data,2)/AccY);


    % load EPI-data 
    EPI_data_raw = meas.data(:,:,:,:,:,:,1:num_shot,:,:,:);

    % load EPI-nav
    EPI_nav = meas.data_phascor1d(:,:,:,:,:,:,1:num_shot,:,:,:);

    % load PAT reference data
    PAT_ref = meas.patrefscan(:,:,:,:,:,:,:,:,:,:);

    % load PAT-nav
    PAT_nav = meas.patrefscan_phascor(:,:,:,:,:,:,:,:,:,:);  


    % Correct the global phase

    EPI_data=zeross(size(EPI_data_raw));
    for ii=1:num_shot
       EPI_data(:,:,:,1,1,1,ii,:,1,:)=msEPI_kyshift_correction(EPI_data_raw(:,:,:,1,1,1,ii,:,1,:), meas.prot, del_ky(ii), -1) ;
    end

    EPI_data_comb=zeross([size(EPI_data,1),size(EPI_data,2),size(EPI_data,3),size(EPI_data,4),size(EPI_data,5),size(EPI_data,6),1,size(EPI_data,8),size(EPI_data,9),size(EPI_data,10)]);


    % Manually shift the ky lines
    for sh = 1:num_shot
       ky_idx = 1 + del_ky(sh) + offset_ky : AccY : s(EPI_data, 2)  + del_ky(sh) + offset_ky;
       ky_idx(ky_idx >=s(EPI_data,2)) = [] 
       EPI_data_comb(:, ky_idx,:,:,:,:,1,:,:,:) =  EPI_data(:,startLine:AccY:startLine+AccY*(length(ky_idx)-1),:,:,:,:,sh,:,:,:);
    end

    % Correct EPI-data
    EPI_data_cor_kspace = ghost_correct_pat_ref_v1_STD_bb(meas.prot, EPI_data_comb, EPI_nav, AccY);


    if diff_dir == 1
        % Correct PAT reference data
        PAT_ref_cor_kspace  = ghost_correct_pat_ref_v1_STD_bb(meas.prot, PAT_ref, PAT_nav, AccY);

        % zero-pad the reference data
        PAT_ref_cor_kspace_zpad = zpad(squeeze(PAT_ref_cor_kspace),nLine,nColumn,nCoil,nSlice);
        PAT_ref_cor_kspace_zpad_ext = zeros(nLine,nColumn,nCoil,1,1,1,1,1,1,nSlice);
        PAT_ref_cor_kspace_zpad_ext(:,:,:,1,1,1,1,1,1,:) = PAT_ref_cor_kspace_zpad;
        PAT_ref_cor_img_ext = mrir_iDFT_freqencode(mrir_iDFT_phasencode(PAT_ref_cor_kspace_zpad_ext));

        PAT_ref_cor_img = single( squeeze(PAT_ref_cor_img_ext) ); 
        PAT_ref_cor_img = PAT_ref_cor_img(1+end/4:3*end/4,:,:,:);

        save PAT_ref_cor_img_subject1 PAT_ref_cor_img 
    end

    EPI_data_cor_img = single(sq(ifft2call(EPI_data_cor_kspace)));
    EPI_data_cor_img = EPI_data_cor_img(1+end/4:3*end/4,:,:,:);

    mosaic(rsos(PAT_ref_cor_img,3),5,8,1,'',[0,.1],90)
    mosaic(rsos(EPI_data_cor_img,3),5,8,2,'',[0,3e-4],90)

    save([file_path, 'EPI_data_cor_img_shots', num2str(shot_accl(1)), 'to', num2str(shot_accl(end)),'.mat'], 'EPI_data_cor_img')
    
end



%--------------------------------------------------------------------------
%% load data
%--------------------------------------------------------------------------

    
addpath /autofs/cluster/kawin/berkin/Matlab_Code/TOOLBOXES/SENSE_LSQR_Toolbox
addpath /autofs/cluster/kawin/berkin/Matlab_Code_New/VC-MUSSELS


recon_path = '/autofs/space/lucifer_001/users/data/2018_11_15_Bay4_InVivo_diffusion_4subjects_NEATR/subject4/';


load([recon_path, 'PAT_ref_cor_img'])
load([recon_path, 'Receive_subject4_R9_9shot'])


% load([recon_path, 'EPI_data_cor_img_shots10to18'])
% load([recon_path, 'EPI_data_cor_img_shots19to27'])
% load([recon_path, 'EPI_data_cor_img_shots28to36'])
% load([recon_path, 'EPI_data_cor_img_shots37to45'])
load([recon_path, 'EPI_data_cor_img_shots46to54'])
% load([recon_path, 'EPI_data_cor_img_shots55to63'])

mosaic( rsos(PAT_ref_cor_img,3),5,8,1,'',[0,.1],90 )
mosaic( rsos(EPI_data_cor_img,3),5,8,2,'',[0,3e-4],90 )


[N(1),N(2),num_chan,num_slc] = size(EPI_data_cor_img)


[meas.prot, meas.evp] = read_meas_prot([recon_path,'meas_MID00333_FID66841_msEPI_1mm_PAT9_9shot.dat']);


%--------------------------------------------------------------------------
%% SENSE & MUSSELS: MB1 
%--------------------------------------------------------------------------


offset_ky = 0;
AccY = meas.prot.lAccelFactPE;           %PAT factor
del_ky = 0 : AccY-1; 


disp_on = 1;
num_ech = 1;
num_shot = 9;

shot_accl = 1:1:num_shot
num_sh = length(shot_accl);


Img_Sense = zeross([N,num_ech,num_sh,num_slc]);
Img_Fista = zeross([N,num_ech,num_sh,num_slc]);
% Receive = zeross([N,num_chan,num_slc]);


for slc_index = 33%1:num_slc
    disp('%--------------------------------------------------------------------------')
    disp(['Slice: ', num2str(slc_index)])
    disp('%--------------------------------------------------------------------------')
    
    
    img_ref = PAT_ref_cor_img(:,:,:,slc_index);                                               % kx, ky, chan, slc
    kspace_slc = sq(fft2call(EPI_data_cor_img(:,:,:,slc_index,1:num_ech)));         % kx, ky, chan, echo, shot, slc
 
    kspace_slice = zeross([N,num_chan,num_ech,num_sh]);

    % shift the ky lines
    for sh = 1:num_sh
       ky_idx = 1 + del_ky(sh) + offset_ky : AccY : N(2)  + del_ky(sh) + offset_ky;
       ky_idx(ky_idx >= N(2)) = []; 

       kspace_slice(:, ky_idx,:,:,sh) = kspace_slc(:,ky_idx,:,:);
    end
    
    %--------------------------------------------------------------------------
    % Espirit
    %--------------------------------------------------------------------------
    
%     c = 0.7;
%     num_acs = 24;
% 
%     receive = sq( CoilSense_ESPIRIT2d( permute( img_ref , [1,2,4,3]), c, min(num_acs), recon_path ) );
% 
%     if disp_on
%         mosaic(rsos(img_ref,3), 1, 1, 1, '', [0,1e-1], 90)    
%         mosaic(coil_combine(img_ref, receive), 1, 1, 2, '', [0,1e-1], 90)    
%     end
    
%     Receive(:,:,:,slc_index) = receive;
    receive = Receive(:,:,:,slc_index);
    
    %--------------------------------------------------------------------------
    % SENSE 
    %--------------------------------------------------------------------------

    lsqr_iter = 200;
    lsqr_tol = 1e-3;

    param = [];
    param.sens = receive;
    param.N = N;
    param.num_chan = num_chan;
    param.lambda = 1e-3;

    img_sense = zeross([N, num_ech, num_sh]);

    for sh = 1:num_sh
        param.m2d = sq(kspace_slice(:,:,:,1,sh)~=0);

        for te = 1:num_ech
            k = kspace_slice(:,:,:,te,sh);

            res = lsqr(@apply_sense_tikc, cat(1, k(:), zeross([prod(N),1])), lsqr_tol, lsqr_iter, [], [], [], param);  

            img_sense(:,:,te,sh) = reshape(res, N);
        end
    end

    Img_Sense(:,:,:,:,slc_index) = img_sense;

    if disp_on
        mosaic(mean(abs(img_sense),4), 1, num_ech, 3, 'sense', [0,4e-4], 90)
    end
    
    
    %--------------------------------------------------------------------------
    % MUSSELS + FISTA
    %--------------------------------------------------------------------------

    fista = 1;

    num_iter = 200;
    tol = 0.2;

    winSize = [7,7];    % local k-space window size  
    lambda = 1.0;       % hard threshold: effective num shots (rank constraint) default=1.0   

    Sens = repmat(receive, [1,1,1,num_sh]);

    cSens = conj(Sens);
    Sens2 = sum(abs(Sens).^2, 3);
    
    img_fista = zeross([N,num_ech,num_sh]);

    for ne = 1:num_ech
        k_slice = sq(kspace_slice(:,:,:,ne,:)); 
        
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

            k_pocs = Row2im(A, [N, num_sh], winSize);
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
            x_kneg1 = x_k;

            if update < tol
                break
            end

            if ~mod(t,10) && disp_on
                mosaic(mean(abs(x_k),4), 1, num_ech, 4, '', [0,4e-4], 90)
            end
        end
        
        img_fista(:,:,ne,:) = x_k;
    end
    
    
    if disp_on
        mosaic(mean(abs(img_fista),4), 1, num_ech, 4, '', [0,4e-4], 90)
    end
    
    Img_Fista(:,:,:,:,slc_index) = img_fista;

    
    if ~mod(slc_index,5)
%         save([recon_path, 'Img_Fista_subject4_R9_9shot_shots55to63.mat'], 'Img_Fista')
%         save([recon_path, 'Img_Sense_subject4_R9_9shot_shots55to63.mat'], 'Img_Sense')
%         save([recon_path, 'Receive_subject4_R9_9shot.mat'], 'Receive')
    end
end

 

%--------------------------------------------------------------------------
%% SENSE & MUSSELS: SMS 
%--------------------------------------------------------------------------

offset_ky = 0;
AccY = meas.prot.lAccelFactPE;           %PAT factor
del_ky = 0 : AccY-1; 

disp_on = 0;
num_shot = 9;

shot_accl = 1:2:num_shot
num_sh = length(shot_accl);


num_ech = 1;
MB_factor = 2;

slc_skip = round(num_slc / MB_factor);


Img_sms_Sense = zeross([N, num_ech, num_sh, num_slc]);
Img_sms_Mussels = zeross([N, num_ech, num_sh, num_slc]);


for slc_select = 15%1:slc_skip
    slc_index = slc_select : slc_skip : num_slc

    disp('%--------------------------------------------------------------------------')
    disp(['Slice: ', num2str(slc_index)])
    disp('%--------------------------------------------------------------------------')

            
    kspace_slc = fft2call( permute( EPI_data_cor_img(:,:,:,slc_index,1:num_ech), [1,2,3,5,4] ) );         % kx, ky, chan, echo, shot, slc
    
    
    kspace_slice = zeross([N, num_chan, num_ech, num_sh, MB_factor]);

    % shift the ky lines
    sh_idx = 1;
    for sh = shot_accl
       ky_idx = 1 + del_ky(sh) + offset_ky : AccY : N(2)  + del_ky(sh) + offset_ky;
       ky_idx(ky_idx >= N(2)) = [];

       kspace_slice(:, ky_idx,:,:,sh_idx,:) = kspace_slc(:,ky_idx,:,:,:);
       sh_idx = sh_idx + 1;
    end
           
    
    % collapsed k-space 
    kspace_colaps = sum(kspace_slice, 6);   
 

    %--------------------------------------------------------------------------
    % Espirit
    %--------------------------------------------------------------------------

%     Receive = zeross([N, num_chan, MB_factor]);

%     c = 0.7;
%     num_acs = 24;
% 
%     for slc = 1:MB_factor
%         im_ref = img_ref(:,:,:,slc);
% 
%         % zero out beginning and end of readout in fleet reference to
%         % avoid supurious mask region outside of brain
%         im_ref(1:10,:,:)=0;
%         im_ref(end-9:end,:,:)=0;
% 
%         receive = sq( CoilSense_ESPIRIT2d( permute( im_ref , [1,2,4,3]), c, min(num_acs), recon_path ) );
% 
%         Receive(:,:,:,slc) = receive;
% 
%         if disp_img
%             mosaic(rsos(im_ref,3), 1, 1, slc, '', genCaxis(rsos(im_ref,3)), 90), setGcf(.5)    
%             mosaic(coil_combine(im_ref, receive), 1, 1, 10+slc, '', genCaxis(rsos(im_ref,3)), 90)    
%         end
%     end
% 
%     Receive = cat(1, Receive(:,:,:,1), Receive(:,:,:,2));
    Rec = cat(1, Receive(:,:,:,slc_index(1)), Receive(:,:,:,slc_index(2)));


    %--------------------------------------------------------------------------
    % slice-SENSE/GRAPPA 
    %--------------------------------------------------------------------------

    lsqr_iter = 200;
    lsqr_tol = 1e-3;

    param = [];
    param.sens = Rec;
    param.N = N .* [2,1];
    param.num_chan = num_chan;
    param.lambda = 1e-3;

    img_sense = zeross([N .* [2,1], num_ech, num_sh]);

    kspc_colaps = kspace_colaps;
    msk = kspc_colaps~=0;                               % Ry in-plane accl

    img_colaps = ifft2call(kspc_colaps);                % collapsed image
    img_colaps = repmat(img_colaps, [2,1,1,1,1]) / 2;   % repmat and scale to emulate Rx undersampling

    kspc_colaps = fft2call(img_colaps);
    kspc_colaps(2:2:end,:,:,:,:) = 0;                   % Rx readout accl
    kspc_colaps = kspc_colaps .* repmat(msk, [2,1,1,1,1]);


    for sh = 1:num_sh
        param.m2d = sq(kspc_colaps(:,:,:,1,sh)~=0);

        for te = 1:num_ech
            k = kspc_colaps(:,:,:,te,sh);

            res = lsqr(@apply_sense_tikc, cat(1, k(:), zeross([prod(N) * MB_factor,1])), lsqr_tol, lsqr_iter, [], [], [], param);  

            img_sense(:,:,te,sh) = reshape(res, N .* [2,1]);
        end
    end

    
    if disp_on
        mosaic(mean(abs(img_sense),4), 1, num_ech, 1, 'sms-sense', [0,4e-4], 90)
    end
    
    Img_sms_Sense(:,:,:,:,slc_index) = cat(5, img_sense(1:end/2,:,:,:), img_sense(1+end/2:end,:,:,:) );

   
    %--------------------------------------------------------------------------
    % SMS-MUSSELS: Fista
    %--------------------------------------------------------------------------

    fista = 1;

    num_iter = 200;
    tol = 0.3;

    winSize = [1,1] * 7;    % local k-space window size  
    lambda = 1.25;             % hard threshold: effective num shots (rank constraint) default=1.0   

    Sens = repmat(Rec, [1,1,1,num_sh]);

    cSens = conj(Sens);
    Sens2 = sum(abs(Sens).^2, 3);

    img_fista = zeross([MB_factor*N(1), N(2), num_ech, num_sh]);

    tic
    for ne = 1:num_ech
        k_slice = sq(kspc_colaps(:,:,:,ne,:)); 

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

            k_pocs = Row2im(A, [MB_factor*N(1), N(2), num_sh], winSize);
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
            x_kneg1 = x_k;

            if update < tol
                break
            end

            if disp_on && ~mod(t,10)
                mosaic(mean(abs(x_k),4), 1, num_ech, 2, 'sms-mussels', [0,4e-4], 90)
            end
        end

        img_fista(:,:,ne,:) = x_k;
    end
    toc
    
    if disp_on
        mosaic(mean(abs(img_fista),4), 1, num_ech, 2, 'sms-mussels', [0,.3], 90)
    end

    Img_sms_Mussels(:,:,:,:,slc_index) = cat(5, img_fista(1:end/2,:,:,:), img_fista(1+end/2:end,:,:,:) );
    
   
    if ~mod(slc_index,2)
%         save([recon_path, 'Img_Fista_subject4_R9_MB2_5shot_shots55to63.mat'], 'Img_sms_Mussels')
%         save([recon_path, 'Img_Sense_subject4_R9_MB2_5shot_shots55to63.mat'], 'Img_sms_Sense')
    end
end
    


%--------------------------------------------------------------------------
%% load data
%--------------------------------------------------------------------------

    
addpath /autofs/cluster/kawin/berkin/Matlab_Code/TOOLBOXES/SENSE_LSQR_Toolbox
addpath /autofs/cluster/kawin/berkin/Matlab_Code_New/VC-MUSSELS


recon_path = '/autofs/space/lucifer_001/users/data/2018_11_15_Bay4_InVivo_diffusion_4subjects_NEATR/subject4/';


load([recon_path, 'PAT_ref_cor_img'])
load([recon_path, 'Receive_subject4_R9_9shot'])


% load([recon_path, 'EPI_data_cor_img_shots10to18'])
% load([recon_path, 'EPI_data_cor_img_shots19to27'])
% load([recon_path, 'EPI_data_cor_img_shots28to36'])
load([recon_path, 'EPI_data_cor_img_shots37to45'])
% load([recon_path, 'EPI_data_cor_img_shots46to54'])
% load([recon_path, 'EPI_data_cor_img_shots55to63'])

mosaic( rsos(PAT_ref_cor_img,3),5,8,1,'',[0,.1],90 )
mosaic( rsos(EPI_data_cor_img,3),5,8,2,'',[0,3e-4],90 )


[N(1),N(2),num_chan,num_slc] = size(EPI_data_cor_img)


[meas.prot, meas.evp] = read_meas_prot([recon_path,'meas_MID00333_FID66841_msEPI_1mm_PAT9_9shot.dat']);



load([recon_path, 'Img_Fista_subject4_R9_9shot_shots37to45'])
load([recon_path, 'Img_Fista_subject4_R9_MB2_5shot_shots37to45'])

% load([recon_path, 'Img_Fista_subject4_R9_9shot_shots46to54'])
% load([recon_path, 'Img_Fista_subject4_R9_MB2_5shot_shots46to54'])

% load([recon_path, 'Img_Fista_subject4_R9_9shot_shots55to63'])
% load([recon_path, 'Img_Fista_subject4_R9_MB2_5shot_shots55to63'])


mosaic(mean(abs(Img_Fista),4),5,8,1,'',[0,4e-4],90)
mosaic(mean(abs(Img_sms_Mussels),4),5,8,2,'',[0,4e-4],90)



%--------------------------------------------------------------------------
%% load recons
%--------------------------------------------------------------------------
  
shot_accl = [1:2:9];
num_sh = length(shot_accl)

[N(1), N(2), num_echo, num_shot, num_slc] = size(Img_Fista);

T1 = mean(abs(Img_sms_Mussels),4);
T2 = mean(abs(Img_Fista(:,:,:,shot_accl,:)),4);

% whole brain rmse:
rmse(T1,T2)     


slc_select = 7;
slc_skip = 20;
slc_index = slc_select : slc_skip : num_slc;



% slice group
msk_receive = Receive(:,:,1,slc_index)~=0;

img_sms_mussels = Img_sms_Mussels(:,:,:,:,slc_index) .* repmat(permute(msk_receive,[1,2,3,5,4]),[1,1,1,num_sh,1]);
img_R1 = Img_Fista(:,:,:,shot_accl,slc_index) .* repmat(permute(msk_receive,[1,2,3,5,4]),[1,1,1,num_sh,1]);


mosaic(mean(abs(img_sms_mussels),4),1,2,1, ['mussels: ', num2str(rmse(mean(abs(img_sms_mussels),4), mean(abs(img_R1),4)))] ,[0,3.5e-4],90)
mosaic(mean(abs(img_R1),4),1,2,2,'',[0,4e-4],90)


%--------------------------------------------------------------------------
%% bm3d + Mussels -> magnitude
%--------------------------------------------------------------------------

addpath /autofs/cluster/kawin/berkin/Matlab_Code_New/NEATR/BM3D/ 


im_bm3d = 0 * mean(img_sms_mussels,4);

sigma_bm3d = 10;

for slc = 1:s(img_sms_mussels,5)
    for te = 1:s(img_sms_mussels,3)
        y = mean(abs(img_sms_mussels(:,:,te,:,slc)),4);

        scl_fctr = max(y(:));

        [~, y_re] = BM3D(1, y / scl_fctr, sigma_bm3d); 

        y_est = y_re * scl_fctr;

        im_bm3d(:,:,te,1,slc) = y_est;
    end
end


mosaic( cat(3, mean(abs(img_sms_mussels(:,:,:,:,1)),4), mean(abs(img_sms_mussels(:,:,:,:,2)),4)), 1, 2, 1, num2str(rmse(t1,t2)), [0,4e-4], 90 )
mosaic( cat(3, mean(abs(im_bm3d(:,:,:,:,1)),4), mean(abs(im_bm3d(:,:,:,:,2)),4)), 1, 2, 2, ['bm3d ', num2str(rmse(im_bm3d,t2))], [0,4e-4], 90 )


% rmse(mean(abs(img_sms_mussels),4),t2)


%--------------------------------------------------------------------------
%% bm3d + Mussels -> complex valued
%--------------------------------------------------------------------------

addpath /autofs/cluster/kawin/berkin/Matlab_Code_New/NEATR/BM3D/ 


img_bm3d = 0 * img_sms_mussels;

sigma_bm3d = 8;

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


img_bm3d = img_bm3d .* repmat(permute(msk_receive,[1,2,3,5,4]),[1,1,1,num_sh,1]);


mosaic( cat(3, mean(abs(img_bm3d(:,:,:,:,1)),4), mean(abs(img_bm3d(:,:,:,:,2)),4)), 1, 2, 3, ['complex bm3d ', num2str(rmse(mean(abs(img_bm3d),4), mean(abs(img_R1),4)))], [0,3.5e-4], 90 )



%--------------------------------------------------------------------------
%% magnitude-valued unet
%--------------------------------------------------------------------------

load /autofs/space/lucifer_001/users/data/2018_11_15_Bay4_InVivo_diffusion_4subjects_NEATR/res_unet_fista_R9_MB2_5shot_patch64x64_100epochs_dwi_5chan
 
Res_unet = permute(Res_unet,[2,4,3,1]);

Res_unet = imrotate(reshape(Res_unet, [N, 1, num_sh, num_slc, s(Res_unet,4)/num_slc]), 180);


dwi_index = 4;


mosaic(mean(abs(Img_Fista),4),5,8,11,'',[0,4e-4],90)
mosaic(mean(abs(Img_sms_Mussels),4),5,8,12,'',[0,4e-4],90)
mosaic(mean(abs(Res_unet(:,:,:,:,:,dwi_index)),4),5,8,13,'',[0,4e-4],90)


res_unet = Res_unet(:,:,:,:,slc_index,dwi_index) .* repmat(permute(msk_receive,[1,2,3,5,4]),[1,1,1,num_sh,1]);


mosaic( cat(3, mean(abs(res_unet(:,:,:,:,1)),4), mean(abs(res_unet(:,:,:,:,2)),4)), 1, 2, 4, ['unet  ', num2str(rmse(mean(abs(res_unet),4), mean(abs(img_R1),4)))], [0,3.5e-4], 90 )



%--------------------------------------------------------------------------
%% scale data
%--------------------------------------------------------------------------

scl_factor = max(max(max(mean(abs(img_R1),4))));

kspc_scale = kspc_colaps / scl_factor;


% % --------------------------------------------------------------------------
% % % phase cycling & jvc-sense
% % --------------------------------------------------------------------------
% % 
% % addpath /cluster/kawin/berkin/Matlab_Code_New/NEATR/
% % addpath /autofs/cluster/kawin/berkin/Matlab_Code_New/Machine_Learning/phase_cycling-master/
% % setPath
% % 
% % 
% % num_ech = 1;
% % use_vc = 1;
% % 
% % phase-cycling param
% % nouteriter = 1;
% % h = 10;
% % ninneriter = 50;
% % 
% % 
% % C = Identity;
% % M = Identity;
% % P = Identity;
% % S = ESPIRiT( Rec );
% % 
% % 
% % ncycles = 16;
% % W = {};
% % 
% % lambda_phase = 1e-3;            % L1 wavelet on phase
% % lambda_phase = 1e-3;            % L1 wavelet on phase
% % Pp = wave_thresh('db4', 3, lambda_phase);
% % 
% % 
% % sense param
% % lsqr_iter = 200;
% % lsqr_tol = 1e-6;
% % lambda_L2 = 1e-2;               % L2 for Sense
% % 
% % cs_outeriter = 2;               % CS outer iters
% % lambda_cs = 1e-6;               % TV for Sense
% % lambda_cs = 0e-6;               % TV for Sense
% % 
% % m2d = sq(kspc_colaps(:,:,:,1,:)~=0);
% % img_jvcsense = zeross([N .* [2,1], num_ech]);
% % 
% % PHS_result = zeross([N .* [2,1], num_ech, num_sh]);
% % 
% % 
% % tic
% % for te = 1:num_ech
% %     magnitude u-net initialization (phase from Mussels)
% %     m = sq( mean(abs( Res_unet(:,:,te,:,slc_index) ), 4) ) / scl_factor;       
% %     x = cat(1, img_sms_mussels(:,:,te,:,1), img_sms_mussels(:,:,te,:,2)) / scl_factor;
% %     
% %     bm3d initialization
% %         m = sq( mean(abs( img_bm3d(:,:,te,:,:) ), 4) );
% %     x = cat(1, img_bm3d(:,:,te,:,1), img_bm3d(:,:,te,:,2));
% %     
% %     mussels initialization
% %     m = sq( mean(abs( Img_sms_Mussels(:,:,te,:,slc_index) ), 4) );
% %     x = cat(1, Img_sms_Mussels(:,:,te,:,slc_index(1)), Img_sms_Mussels(:,:,te,:,slc_index(2)));
% %         
% % 
% %     m = cat(1, m(:,:,1), m(:,:,2)) .* (Rec(:,:,1)~=0);
% %     
% %     Phs = sq( angle( x ) ) .* repmat( Rec(:,:,1)~=0, [1,1,num_sh] );
% %     Phs_result = 0 * Phs;
% %     
% %     for n = 1:nouteriter
% %         Mm = M * m;
% %         alphap = 1.0 / lipschitz(P) / (max(abs(Mm(:)))^2 + eps) * h;
% % 
% %         for t = 1:num_sh
% %             y = sq( kspc_colaps(:,:,:,te,t) ) / scl_factor;
% %             F = p2DFT(m2d(:,:,:,t), [N .* [2,1], num_chan]);                 
% %                  
% %             initial guess: from Mussels recon
% %             if n > 1
% %                 p0 = sq(PHS_result(:,:,te,t));
% %             else
% %                 p0 = Phs(:,:,t);
% %             end
% %             p = p0;          
% %                         
% %             for c = 1:ncycles
% %                 p_c = angle( exp(1i * p0) .* exp(1i * (c-1) * 2 * pi /  ncycles) );
% %                 
% %                 if c == 1
% %                     p0_c = p_c;
% %                 end
% %                 W{c} = p_c - p0_c;
% %             end
% %             
% %             for itp = 1:ninneriter
% %                 if ~mod(itp,50)
% %                     disp(num2str(itp))
% %                 end
% %                 
% %                 if isempty(W)
% %                     w = 0;
% %                 else
% %                     ti = randi(length(W));
% %                     w = W{ti};
% %                 end
% %                 
% %                 expPp = exp(1j * (P * p));
% %                 r = C' * (S' * (F' * (y - (F * (S * (C * (Mm .* expPp)))))));
% %                 
% %                 if lambda_phase
% %                     p = Pp(p + alphap * real(-1j * (P' * (conj(Mm) .* conj(expPp) .* r))) + w, alphap) - w;                    
% %                 else
% %                     p = p + alphap * real(-1j * (P' * (conj(Mm) .* conj(expPp) .* r)));
% %                 end        
% %             end
% %             
% %             p = p .* (Rec(:,:,1)~=0);
% % 
% %             Phs_result(:,:,t) = p;            
% %             PHS_result(:,:,te,t) = p;
% % 
% %             mosaic( cat(2, cat(1, p0(1:end/2,:), p(1:end/2,:), 1*(p(1:end/2,:)-p0(1:end/2,:))), cat(1, p0(1+end/2:end,:), p(1+end/2:end,:), 1*(p(1+end/2:end,:)-p0(1+end/2:end,:)))),...
% %                 1, 1, 90*nouteriter+t, 'p0, pEst, 5xdiff', [-pi,pi], 90), setGcf(.2)
% %             
% %             pause
% %         end
% % 
% %         
% %         enforce real valued recon using conjugate symmetric data
% %         kspace_joint = zeross([N .* [2,1], num_chan * num_sh]);
% %         M2d = zeross([N .* [2,1], num_chan * num_sh]);
% %         sens = zeross([N .* [2,1], num_chan * num_sh]);
% % 
% %         for t = 1:num_sh
% %             kspace_coils = sq( kspc_colaps(:,:,:,te,t) ) / scl_factor;
% %             kspace_joint(:,:, 1 + (t-1) * num_chan : t * num_chan ) = kspace_coils;
% % 
% %             M2d(:,:, 1 + (t-1) * num_chan : t * num_chan ) = m2d(:,:,:,t);
% %             sens(:,:, 1 + (t-1) * num_chan : t * num_chan ) = dot_mult( Rec, exp(1i * Phs_result(:,:,t)) );
% %         end
% % 
% %         if use_vc
% %             virtual coil data
% %             temp = fft2c(conj(ifft2c(kspace_joint~=0))) > 1e-6;
% %         
% %             Kspace_joint = fft2c(conj(ifft2c(kspace_joint))) .* temp;
% %             Kspace_joint = cat(3, kspace_joint, Kspace_joint);           
% %             
% %             Rec_use = cat(3, sens, conj(sens));
% %         else
% %             Kspace_joint = kspace_joint;
% %             Rec_use = sens;
% %         end
% %         
% %         if lambda_cs > 0
% %             CS for JVC-Sense
% %             param = init_cs;
% %             param.data = Kspace_joint;
% % 
% %             param.TV = TV2D;
% %             param.M3d = Kspace_joint~=0;
% % 
% %             param.TVWeight = lambda_cs;
% %             param.L1Weight = 0*lambda_cs;
% %             
% %             param.num_chan = s(Kspace_joint,3);
% %             param.Receive = Rec_use;
% %             param.Ct = conj(param.Receive);
% % 
% %             m = zeross(N .* [2,1]);
% % 
% %             for n_cs = 1:cs_outeriter
% %                 m = fnlCg_l1_sense(m, param);
% %             end            
% %         else
% %             param = [];
% %             param.sens = Rec_use;
% % 
% %             param.N = N .* [2,1];
% %             param.num_chan = size(param.sens,3);
% % 
% %             param.m2d = Kspace_joint~=0;
% %             param.lambda = lambda_L2;
% % 
% %             res = lsqr(@apply_sense_tikc, cat(1, Kspace_joint(:), zeross([prod(N) * 2, 1])), lsqr_tol, lsqr_iter, [], [], [], param);  
% %             
% %             m = abs( reshape(res, N.*[2,1]) );
% %         end
% %     end
% %         
% %     img_jvcsense(:,:,te) = m; 
% % end
% % toc
% % 
% % 
% % im_jvcsense = cat(3, img_jvcsense(1:end/2,:), img_jvcsense(1+end/2:end,:)) * scl_factor;
% % 
% % 
% % mosaic( im_jvcsense, 1, 2, 6, ['jvc sense  ', num2str(rmse(im_jvcsense, mean(abs(img_R1),4), 1))], [0,3.5e-4], 90 )

 

