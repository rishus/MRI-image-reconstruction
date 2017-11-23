
% sets up the following variables for the calling program:
% 1. reconstruction points r and their counterparts t.
% 2. quadrature points F.u for formula 18 (each t leads to a different set of u's).
% 3. what order kernel do we want?: maxorder
% 4. the filter: F.H 


tol = 10^(-10); % tolerance for quadrature interval length, num eval of K etc.
maxorder = 4; R = 1.0; eps = 0.001;
Nr = 2^5; Nfrier = 2^4;
F = struct('Nr', Nr, 'Nfrier',Nfrier);


% requisite reconstruction points:
% (i) r values
r = linspace(0.01,0.9,Nr);
% (ii) phi values
Nphi = F.Nfrier;
F.phi = (2*pi/Nphi)*(0:Nphi-1)';  % leave out 2*pi


% for each t = R-r, we have a bunch of u values (quad. pts.)
% these will be needed for G_n(u)'s later in the main program
F.u = cell(Nr,1);  
F.H = cell(Nr,1);  % that's the goal here: calculate and store H_n(t,u)

t = R - r;
maxNu = 0;
for it =  0:F.Nr-1
  [gquad,wts] = gauss_quad_rules(max(4,16-floor(it/2)));
  arg = 0.5*(t(it+1)-2.0*eps)*gquad + 0.5*(t(it+1)+2.0*eps);
  F.u{it+1,3} = length(arg);  % no. of quadrature points used for this t.
  Nu = F.u{it+1,3};
  if (Nu > maxNu)
	maxNu = Nu
  endif
end
 
	
Lhandle = zeros(maxNu,maxorder,Nfrier,Nr);
for it =  0:F.Nr-1
  printf("For iteration = %d \n", it);
  [gquad,wts] = gauss_quad_rules(max(4,16-floor(it/2)));
  arg = 0.5*(t(it+1)-2.0*eps)*gquad + 0.5*(t(it+1)+2.0*eps);
  F.u{it+1,1} = arg;
  F.u{it+1,2} = wts';  % wts to be eventually used in equation 18 integral.
  F.u{it+1,3} = length(arg);  % no. of quadrature points used for this t.
  Nu = F.u{it+1,3};
  for ifrier = 0:F.Nfrier-1
    for iu = 0:Nu-1
      for order = 0:maxorder-1  % due to C++ communication
	Lhandle(iu+1,order+1,ifrier+1,it+1) = add_arguments([F.u{it+1,1}(iu+1), t(it+1), ifrier, R, order]);
      end
    end
  end
end
  %%%%%%%%%%%%%%%%%%%%%%
  compute_all_integrations();
  %%%%%%%%%%%%%%%%%%%%%%%
for it =  0:F.Nr-1
  Nu = F.u{it+1,3};
  for ifrier = 0:F.Nfrier-1
    for iu = 0:Nu-1
      for order = 0:maxorder-1  % due to C++ communication
	Lres(iu+1,order+1,ifrier+1,it+1) = get_res([Lhandle(iu+1,  order+1, ifrier+1,it+1)]);
      end
    end
  end
end

% for it = 0:F.Nr-1
%   Nu = F.u{it+1,3};
%   for ifrier = 0:F.Nfrier-1
%      L= zeros(Nu,maxorder);
%      for iu = 0:Nu-1
%         for order = 0:maxorder-1  
%            L(iu+1,order+1) = get_res([Lhandle(iu+1,  order+1, ifrier+1,it+1)]);
%         end
%      end
     
%      % update H: add all the columns in L
%      Htmp = zeros(Nu,1);
%      for order = 0:maxorder-1
%          Htmp = Htmp + ((-1)^(order+1))*L(:,(order+1));
%      end
%      F.L{it+1}(:,:,ifrier+1) = L;  % each col cors to a fixed order, L_0, L_1, etc.
%      F.H{it+1}(:,ifrier+1) = Htmp(:);
%      clear Htmp L
%   end
% end

clear Nr Nfrier Nu count dt du order count 
clear iu it ifrier h_wait arg gquad wts
save kernel_4thorder.mat

