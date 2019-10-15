function [D,X] = k_svd(img, D, K, T, eps, iter, C)
%img : Í¼Ïñpatch¾ØÕó
%D  : dictionary of atoms of size n ¡Á K
%T  : sparsity threshold
%eps: error tolerance
%K  : number of atoms in image layer dictionary
%N  : number of k-svd iterations to run
%X  : Sparse representation of Y
 tic;
 n = size(D, 1); %patchÊı
 M = size(img, 2); %imgÁĞÊı
 X = zeros(K, M); 
 for i=1:iter
     if (nargin<7)
         for l=1:M
             X(:,l)=OMP(img(:,l),D,T,eps);
         end
     else
         X = C;
     end
     for j=1:K
         idx = find(X(j, :));
         if isempty(idx)
             E=img - D * X;
			 check = sum(E.^2, 1);
			 ind = find(check == max(check),1);
			 X(j, :)=zeros(1, M);
			 D(:, j)=img(:, ind)/(norm(img(:, ind) + 0.001));
         else
             X(j, :)=zeros(1, M);
			 E= img - D * X;
			 Ek = E(:, idx);
			 [U, S, V]=svd(Ek);
			 X(j, idx)=S(1, 1)*V(:, 1)'; 
			 D(:, j)=U(:, 1)/(norm(U(:, 1) + 0.001));
         end
     end
     toc;
 end
end

