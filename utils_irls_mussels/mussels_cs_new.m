function [x] = mussels_cs_new(b,Sen,SM,ksize,Niter,iter,Lambda,r, tol, x_ini, phase_update, vcc_flag)

% Based on original mussels_cs.m

%   INPUTS:
%   b: k-space data of size Nx x Npf x Nch X Nsh
%   Sen: coil sensitivity maps of size Nx x Ny x Nch x MB
%   SM : Sampling pattern of the shots Nx x Ny x Nch x MB
%   iter : number of outer iterations
%   Niter  : number of inner iterations of pcg
%   ksize: filter size
%   rank: the threshold of hard thresholding
%   Eg: tic;recon= mussels_cs(b,Sen,SM,[6,6],5,5,5e4,1,0);toc

if ~exist('vcc_flag', 'var')
    vcc_flag=1;
end

[N1,N2,~,Nsh]=size(b);
hanwin0=single(window2(N1, N2, @hann));

%%% define the Hankel operator and its inverse
T2 = @ (z) getHankel(z,ksize);
T2t= @ (z) invHankel(z,N1,N2,Nsh,ksize);
Nsz=ksize(1)*ksize(2)*Nsh;


C = @ (z,tflag) afun_pcg(z,Sen,SM,tflag);
x0=C(b,'transp'); clear b

if exist('x_ini', 'var') && ~isempty(x_ini)
    x=reshape(x_ini, [N1,N2,Nsh]);
else
    x=reshape(x0,N1,N2,Nsh);
end
if ~exist('tol', 'var') || isempty(tol)
    tol=1e-10;
end

if vcc_flag
    for int=1:Nsh
        x4(:,:,int)=conj(rot90(squeeze(x(:,:,int)),2));
    end
    D=cat(2,T2(x),T2(x4));
else
    x4=[];
    D=T2(x);
end

[U,S,~] = svd(D'*D,'econ');

if exist('x_ini', 'var') && ~isempty(x_ini)
    sub_img=ifft2c(x_ini);
    phasemap=sub_img./(abs(sub_img)+eps);
    eps0=0;
    sub_img_old=mean(abs(sub_img), 3);
else
    eps0=1*S(r+1,r+1)^0.5;
    phasemap=ifft2c(x.*hanwin0);
    phasemap=phasemap./(abs(phasemap)+eps);
    sub_img_old=abs(mean(ifft2c(x).*conj(phasemap), 3));
end

for in =1: iter
    C = @ (z,tflag) afun_pcg(z,Sen,SM,tflag);
    
    % ================================
    % Nuclear norm implemented as IRLS
    % --------------------
    %  ||Ax-b||^2 + Lambda* ||T(x)||*
    % =||Ax-b||^2 + Lambda* ||T(x)Wn^0.5||^2_F
    % W=(T2'T2)^-0.5
    % --------------------

    s=diag(S+eps0).^(-0.5);
    W=(U*diag(s));
   
    % ================================
    %  PCG: CG Sense update
    % 2A'Ax+2Lambda*T2'(T2(x)Wn^0.5Wn'^0.5)=2A'b
    % 2A'Ax+2Lambda*G'(W)G(W)x=2A'b
    % 2A'Ax+2lambda*deconv(conv(f,x),f')=2A'b
    % 2A'Ax+2lambda*IFFT(FFT(IFFT(FFT(f')*FFT(x'))),FFT(f'))=2A'b
    % ================================
    
    W=W*W';
    
    WWt = @ (z) getWWt(z,W,T2,T2t,Nsz, vcc_flag);
    
    H1=@(x)C(x,'notransp')+(Lambda*reshape(WWt(reshape((x),N1,N2,Nsh)),[],1));
    
    
    [x1,~,relres,Nend,resvec] = pcg(H1,2.*x0(:),1e-10,Niter,[],[],x(:));
    x=reshape((x1),N1,N2,Nsh);

    if exist('x_ini', 'var') && ~isempty(x_ini)
        if exist('phase_update', 'var') && phase_update>0
            phasemap=ifft2c(x.*hanwin0);
            phasemap=phasemap./(abs(phasemap)+eps);
        end
        sub_img=abs(mean(ifft2c(x).*conj(phasemap), 3));
       
        relres1=norm(sub_img(:)-sub_img_old(:))/norm(sub_img_old(:));
        sub_img_old=sub_img;
        if relres1 <= tol
            break;
        end
        
        x=fft2c(sub_img.*phasemap);
    else
        phasemap=ifft2c(x.*hanwin0);
        phasemap=phasemap./(abs(phasemap)+eps);
        sub_img=abs(mean(ifft2c(x).*conj(phasemap), 3));
        
        relres1=norm(sub_img(:)-sub_img_old(:))/norm(sub_img_old(:));
        sub_img_old=sub_img;
        if relres1 <= tol
            break;
        end
        
        x=fft2c(sub_img.*phasemap);
    end

    relres_all(in)=relres1;
    
    if vcc_flag
        for int=1:Nsh
            x4(:,:,int)=conj(rot90(squeeze(x(:,:,int)),2));
        end
        D=cat(2,T2(x),T2(x4));
    else
        x4=[];
        D=T2(x);
    end
    
    
    [U,S,~] = svd(D'*D,'econ');
end
fprintf('IRLS MUSSELS converged at iteration %d\n', in);
end

%------------------- pcg ----------------%
function y = afun_pcg(data,Sen,SM,transp_flag)
CSen=conj(Sen);
N=size(Sen);

if strcmp(transp_flag,'transp')
    
    for j=1:size(data,4) % Nint
        for i=1:N(3) % Nchannel
            idata_tmp1(:,:,i)=sqrt((N(1)*N(2))).*CSen(:,:,i).*fftshift(ifft2(fftshift(data(:,:,i,j))));%.*mask(:,:,j)
        end
        y(:,:,j)=fftshift(fft2(fftshift(sum(idata_tmp1(:,:,:),3))))./sqrt((N(1)*N(2)));
    end
    
elseif strcmp(transp_flag,'notransp')
    
    data=reshape(data,N(1),N(2),[]);
    
    for j=1:size(data,3) % Nint
        gdata = sqrt((N(1)*N(2))).*fftshift(ifft2(fftshift(data(:,:,j))));
        for i= 1:size(Sen,3)% Nchannel
            gdata1=fftshift(fft2(fftshift(gdata.*Sen(:,:,i)))).*SM(:,:,j)./sqrt((N(1)*N(2)));%
            gdata2(:,:,i)=sqrt((N(1)*N(2))).*fftshift(ifft2(fftshift(gdata1))).*CSen(:,:,i);
        end
        gdata4(:,:,j)=fftshift(fft2(fftshift(sum(gdata2,3))))./sqrt((N(1)*N(2)));
    end
    y=2.*gdata4(:);
    
end
end



function fwd_Im = getWWt(x,W,T2,T2t,Nsz, vcc_flag)

if vcc_flag
    for int=1:size(x,3)
        x4(:,:,int)=conj(rot90(squeeze(x(:,:,int)),2));
    end
    D=cat(2,T2(x),T2(x4));
else
    x4=[];
    D=T2(x);
end

D=D*W;
x1=T2t(D(:,1:Nsz));

if vcc_flag
    x2=T2t(D(:,Nsz+1:end));
    for int=1:size(x,3)
        fwd_Im(:,:,int)=x1(:,:,int)+conj(rot90(squeeze(x2(:,:,int)),2));
    end
else
    for int=1:size(x,3)
        fwd_Im(:,:,int)=x1(:,:,int);
    end
end

end

%-----------------------------------%
