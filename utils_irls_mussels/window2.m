% This function creates a 2 dimentional window for a sample image, it takes
% the dimension of the window and applies the 1D window function
% This is does NOT using a rotational symmetric method to generate a 2 window
%
% Disi A ---- May,16, 2013
%     [N,M]=size(imgage);
% ---------------------------------------------------------------------
%     w_type is defined by the following 
%     @bartlett       - Bartlett window.
%     @barthannwin    - Modified Bartlett-Hanning window. 
%     @blackman       - Blackman window.
%     @blackmanharris - Minimum 4-term Blackman-Harris window.
%     @bohmanwin      - Bohman window.
%     @chebwin        - Chebyshev window.
%     @flattopwin     - Flat Top window.
%     @gausswin       - Gaussian window.
%     @hamming        - Hamming window.
%     @hann           - Hann window.
%     @kaiser         - Kaiser window.
%     @nuttallwin     - Nuttall defined minimum 4-term Blackman-Harris window.
%     @parzenwin      - Parzen (de la Valle-Poussin) window.
%     @rectwin        - Rectangular window.
%     @taylorwin      - Taylor window.
%     @tukeywin       - Tukey window.
%     @triang         - Triangular window.
%
%   Example: 
%   To compute windowed 2D fFT
%   [r,c]=size(img);
%   w=window2(r,c,@hamming);
% 	fft2(img.*w);

function w=window2(N,M,w_func, L)

wc=window(w_func,N);
wr=window(w_func,M);
[maskr,maskc]=meshgrid(wr,wc);
w=maskr.*maskc;

%% daiep 161007
% if nargin < 4
%     L=16;
% end
% win0=window(w_func,L);
% if L<N
%     wc=ones(N,1);
%     wc(1:L/2)=win0(1:L/2);
%     wc(end-L/2+1:end)=win0(L/2+1:L);
% else
%     wc=window(w_func,N);
% end
% 
% if L<M
%     wr=ones(M,1);
%     wr(1:L/2)=win0(1:L/2);
%     wr(end-L/2+1:end)=win0(L/2+1:L);
% else
%     wr=window(w_func,M);
% end
% [maskr,maskc]=meshgrid(wr,wc);
% w=maskr.*maskc;

end

