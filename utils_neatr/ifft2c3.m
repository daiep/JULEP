function res = ifft2c3(x)
fctr = size(x,1)*size(x,2);

res = zeros(size(x));


for k = 1:size(x,5)
    for m = 1:size(x,4)
        for n = 1:size(x,3)
            res(:,:,n,m,k) = sqrt(fctr)*fftshift(ifft2(ifftshift(x(:,:,n,m,k))));
        end
    end
end