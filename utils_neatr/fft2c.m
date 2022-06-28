function im = fft2c(d)
% Function performs a centered fft2

im = 0 * d;

scl = sqrt(numel(d(:,:,1,1)));

for c = 1:size(d,4)
    for t = 1:size(d,3)
        im(:,:,t,c) = fftshift(fft2(fftshift(d(:,:,t,c)))) / scl;
    end
end