function cernel = compressKernel(kernel)

% cernel = compressKernel(kernel)
%
% The function takes a kernel generated by dat2Kernel or by dat2Ata and reduces
% the number of kernels by performing inner products. The size of each kernel is
% size 2x-1 where x is the kernel window size before. Effectively this operation
% is equivalent to performing G'G in the image domain. For more information
% 
% (c) Michael Lustig 2010


cernel = zeros(size(kernel,1)*2-1,size(kernel,2)*2-1,size(kernel,3),size(kernel,3));
for n=1:size(kernel,3)
    tmp1 = kernel(:,:,n,:);
    for m=1:size(kernel,3);
        tmp2 = kernel(:,:,m,:);
        tmp3 = zeros(size(kernel,1)*2-1,size(kernel,2)*2-1,1,size(kernel,4));
        for l=1:size(kernel,4)
            tmp3(:,:,1,l) = filter2(conj(tmp1(:,:,1,l)),(tmp2(:,:,1,l)),'full');
        end
        cernel(:,:,n,m) = sum(tmp3,4);
        
    end
end
