load  kernel_4thorder.mat

tol = 10^(-10); % tolerance for quadrature interval length, num eval of K etc.
maxorder = 4; R = 1.0; eps = 0.001;
Nr = 2^5; Nfrier = 2^4;
r = linspace(0.01,0.9,Nr);
Nphi = F.Nfrier;

Nphi
F.Nr

for it = 0:F.Nr-1
  Nu = F.u{it+1,3};
  for ifrier = 0:F.Nfrier-1
     L= zeros(Nu,maxorder);
     for iu = 0:Nu-1
        for order = 0:maxorder-1  
           L(iu+1,order+1) = Lres(iu+1,  order+1, ifrier+1,it+1);
           Lres(iu+1,  order+1, ifrier+1,it+1);
        end
     end
     
     % update H: add all the columns in L
     Htmp = zeros(Nu,1);
     for order = 0:maxorder-1
         Htmp = Htmp + ((-1)^(order+1))*L(:,(order+1));
     end
     F.L{it+1}(:,:,ifrier+1) = L;  % each col cors to a fixed order, L_0, L_1, etc.
     F.H{it+1}(:,ifrier+1) = Htmp(:);
     clear Htmp L
  end
end
F.H

