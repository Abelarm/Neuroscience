%[Ealpha, Ebeta, Omega, b, Evarphi, tau, nu] = vb_classification_dirichlet(x,t,hyper,Kdef,bias,display,MAXIT,KATol,SampPars,Ksupplied)
%
%Heterogenous Kernel Combination - Variational Bayes Classification
%Copyright: Simon Rogers & Mark Girolami, March 2004
%
%INPUTS
%x - Input data points (columns = features)
%t - Targets (+-1)
%hyper - Structure holding hyper-parameter values
%Kdef - Structure holding kernel definitions
%   Kdef(1).ktype = 'supplied' if Kernel is supplied by user.([N M K]
%   array)
%   Otherwise:  Kdef(i).ktype = {'rbf','linear','polynomial'}
%               Kdef(i).kpar = sigma for rbf, order for polynomial
%bias: 1 for bias, 0 otherwise (ignored if kernel is supplied by user)
%display: 0 for no feedback, 1 for iteration number, 2 for plots
%MAXIT: Maximum number of iterations
%KATol: Tolerance for change in Kernel*Alpha
%SampPars - Structure of sampling parameters
%   SampPars.NoS: Number of samples
%   SampPars.MAXIT: Number of sampling iterations
%   SampPars.TOL: Sampling tolerance
%Ksupplied: supplied kernel (optional)
%
%
%
%OUTPUTS
%Ealpha - Expected value of \alpha
%Ebeta - Expected value of \beta
%Omega,b,Evarphi,tau,nu - Optional .. see paper for details
%
%
%
function [Ealpha, Ebeta, Omega, b, Evarphi, tau, nu] = vb_classification_dirichlet(x,t,hyper,Kdef,bias,display,MAXIT,KATol,SampPars,Ksupplied)
if display > 0
  fprintf('\nHeterogenous Kernel Combination\nVB Classification\n--------------------------');
  fprintf('\nCreating Kernel: ');
end


%Create the kernel is it has not been supplied
switch Kdef(1).ktype
    case 'supplied'
        Kall = Ksupplied;
        [N,M,K] = size(Kall);
    otherwise
        [N,M] = size(x);
        M = N;
        if bias > 0
            M = M + 1;
        end
        K = length(Kdef);
        Kall = zeros(N,M,K);
        for k = 1:K
           temp = kernel(x,x,Kdef(k).ktype,Kdef(k).kpar); 
            q = (1./diag(temp));
            temp = temp.*repmat(q',N,1);
            if bias > 0
                temp = [temp repmat(1,N,1)];
            end
            Kall(:,:,k) = temp;
        end
end

%Hyper-paramaters are stored in hyper
%eg. sigma = hyper.sigma;

%Set-up Parameters
%Alpha mean and covariance
Ealpha = zeros(M,1);
Salpha = zeros(M);
Ealpha_alphaT = eye(M);
%Phi moments
Ephi = repmat(hyper.sigma/hyper.varsigma,M,1);
Elog_phi = zeros(M,1);
%Beta
Ebeta = repmat(1/K,K,1);
Ebeta_betaT = repmat(1,K,K);
Elog_beta = log(Ebeta);
%Var
Evarphi = repmat(1,K,1);
Elog_varphi = log(Evarphi);
Elog_gamma_varphi = log(gamma(Evarphi));
Elog_gamma_sum_varphi = log(gamma(sum(Evarphi)));
%Xi
Xi = repmat(1,N,1);


%Make current matrix
thisK = zeros(N,M);
for k = 1:K
    thisK = thisK + Ebeta(k)*Kall(:,:,k);
end

%Initial Xi value
Xi = sqrt(diag(thisK*Ealpha_alphaT*thisK'));
Xi(find(Xi<1e-10)) = 1e-10;
Nab = diag(2*Lam(Xi));

%Diagnostics
Gall = [];%Evolution of Gamma
KAChange = [];KAOld = 1;%Evolution of Output function (i.e. thisK * Alpha)
Lall = [];
Aall = [];
boundall = [];
%Subplot sizes for display
maxpr = 4;maxpc = 2;

if display > 0
  fprintf('\nIteration: ');
end

for it = 1:MAXIT
    if display > 0
        temp = it;
        while temp > 1
            fprintf('\b');
            temp = temp/10;
        end
        fprintf('%g',it);
    end

   ppos = 1;  
   
   TempA = zeros(M);
   for i = 1:K
       for j = 1:K
            TempA = TempA + Ebeta_betaT(i,j)*(Kall(:,:,i)'*Nab*Kall(:,:,j));
       end
   end
   TempA = TempA + diag(Ephi);
   Salpha = inv(TempA);
   Ealpha = (0.5) * Salpha * (thisK' * t);
   Ealpha_alphaT = Salpha + Ealpha * Ealpha'; 
   KAChange = [KAChange;norm(thisK*Ealpha - KAOld)];
   KAOld = thisK*Ealpha;
  
   Xi = sqrt(diag(thisK*Ealpha_alphaT*thisK'));
   Xi(find(Xi<1e-10)) = 1e-10;
   Nab = diag(2*Lam(Xi));

   %Omega - needed often
   Omega = zeros(K);
   for i = 1:K
       for j = 1:K
        Omega(i,j) = trace(Nab.*(Kall(:,:,i)*Ealpha_alphaT*Kall(:,:,j)'));
       end
   end

   %Update Phi
   Ephi = (1+2*hyper.sigma)./(diag(Ealpha_alphaT)+2*hyper.varsigma);

   
   if display>1
       subplot(maxpr,maxpc,ppos),bar(Ealpha);ppos = ppos + 1;
       subplot(maxpr,maxpc,ppos),bar(Ephi);ppos = ppos + 1;
       drawnow;
   end
   
   %Update varphi statistics with sampling
   Zphi = 0;%Estimate of normalisation factor
   SumW = 0;
   Offset = 0;


   Evarphi = repmat(0,K,1);
   Elog_varphi = repmat(0,K,1);
   Elog_gamma_varphi = repmat(0,K,1);
   Elog_gamma_sum_varphi = 0;
   Evarphi_old = Evarphi;
   for s_it = 1:SampPars.MAXIT
       [Evarphi, Elog_varphi, Elog_gamma_varphi, Elog_gamma_sum_varphi, SumW, Offset] = sample_varphi(...
           SampPars.NoS, Evarphi, Elog_varphi, Elog_gamma_varphi, Elog_gamma_sum_varphi,...
           SumW,hyper.tau,hyper.nu,Elog_beta,Offset);
       if norm(Evarphi - Evarphi_old)<SampPars.TOL
           break
       end
       Evarphi_old = Evarphi;
   end
   Zphi = SumW/(s_it*SampPars.NoS); 

   
   %Update Beta with sampling
   
   b = repmat(0,K,1);
   for k = 1:K
       b(k) = sum(t.*(Kall(:,:,k)*Ealpha));
   end

   if K>1
     Zbeta = 0;
     Offset = 0;
     SumW = 0;
     Ebeta = repmat(0,K,1);
     Elog_beta = repmat(0,K,1);
     Ebeta_betaT = repmat(0,K,K);
     Ebeta_old = Ebeta;
   
     for s_it = 1:SampPars.MAXIT
       [Ebeta, Elog_beta, Ebeta_betaT, SumW, Offset, eflag] = sample_beta(...
           SampPars.NoS, Ebeta, Elog_beta, Ebeta_betaT, SumW, Evarphi, Omega, b, Offset);
       if norm(Ebeta - Ebeta_old)<SampPars.TOL
         break
       end
       if eflag == -1
         break
       end
       Ebeta_old = Ebeta;
     end
   end
   %Beta contribution to the bound (Normalising Constant)
   Zbeta = SumW/(s_it*SampPars.NoS);

   
 
   
   %Make thisK
   thisK = zeros(N,M);
   for k = 1:K
       thisK = thisK + Ebeta(k)*Kall(:,:,k);
   end
   %Calculate Z
   Z = zeros(K,N);
   for k = 1:K
       Z(k,:) = (Kall(:,:,k)*Ealpha)';
   end
   
   %Hyperparameter updates
   %Sigma
   if hyper.update == 1
        top_const = -log((1/M)*sum(Ephi))+(1/M)*sum(Elog_phi);
        hyper.sigma = update_hyper(hyper.sigma,top_const);
        hyper.varsigma = (M*hyper.sigma)/(sum(Ephi));
        %Tau
        if K>1
          top_const = -log((1/K)*sum(Evarphi))+(1/K)*sum(Elog_varphi);
          hyper.tau = update_hyper(hyper.tau,top_const);
          hyper.nu = (K*hyper.tau)/(sum(Evarphi));
        end
   end
   bound = compute_bound(Ephi,Elog_phi,Ealpha,Ealpha_alphaT,Salpha,...
        hyper,Evarphi,Elog_varphi,Zphi,Ebeta,Elog_beta,Zbeta,...
        Ebeta_betaT,Elog_gamma_varphi,N,K,Elog_gamma_sum_varphi,Omega,b,Xi,Z,t);
   bound = bound + Offset;
   boundall = [boundall;bound];
   if display>1
    subplot(maxpr,maxpc,ppos),semilogy(boundall);ppos = ppos + 1;
    subplot(maxpr,maxpc,ppos),bar(Ebeta);ppos = ppos + 1;title('Beta');
    subplot(maxpr,maxpc,ppos),bar(Evarphi);ppos = ppos + 1;title('VarPhi');
  end
end


ppos = 1;  

TempA = zeros(M);
for i = 1:K
  for j = 1:K
    TempA = TempA + Ebeta_betaT(i,j)*(Kall(:,:,i)'*Nab*Kall(:,:,j));
  end
end
TempA = TempA + diag(Ephi);
Salpha = inv(TempA);
Ealpha = (0.5) * Salpha * (thisK' * t);
Ealpha_alphaT = Salpha + Ealpha * Ealpha'; 
KAChange = [KAChange;norm(thisK*Ealpha - KAOld)];
KAOld = thisK*Ealpha;

Xi = sqrt(diag(thisK*Ealpha_alphaT*thisK'));
Xi(find(Xi<1e-10)) = 1e-10;
Nab = diag(2*Lam(Xi));
%Omega - needed often
Omega = zeros(K);
for i = 1:K
  for j = 1:K
    Omega(i,j) = trace(Nab.*(Kall(:,:,i)*Ealpha_alphaT*Kall(:,:,j)'));
  end
end
b = repmat(0,K,1);
for k = 1:K
  b(k) = sum(t.*(Kall(:,:,k)*Ealpha));
end

tau = hyper.tau;
nu = hyper.nu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = update_hyper(h,top_const)
hOld = h;
tol = 1e-4;
for it = 1:100
    h = hOld*exp((log(hOld) - psi(0,hOld) + top_const)/(-1+hOld*psi(1,hOld)));
    if norm(hOld-h)<tol
        break
    end
    if h == 0
      h = 1;
      break
    end
    if isnan(h)
      h=1;
      break;
    end
    if isinf(h)
      h=1;
      break;
    end
    hOld = h;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bound = compute_bound(Ephi,Elog_phi,Ealpha,Ealpha_alphaT,Salpha,...
    hyper,Evarphi,Elog_varphi,Zphi,Ebeta,Elog_beta,Zbeta,...
    Ebeta_betaT,Elog_gamma_varphi,N,K,Elog_gamma_sum_varphi,Omega,b,Xi,Z,t)
%Alpha
bound = 0.5*sum(Elog_phi) - 0.5*sum(Ephi.*diag(Ealpha_alphaT));
bound = bound - 0.5*safe_log(det(Salpha));
%Phi
bound = bound -(N+1)*log(hyper.sigma) + (N+1)*hyper.sigma*log(hyper.varsigma) + ...
       (1-hyper.sigma)*sum(Elog_phi) - hyper.varsigma * sum(Ephi);
bound = bound - ((hyper.sigma + 0.5)*sum(log(hyper.varsigma + diag(Ealpha_alphaT))) + ...
       (hyper.sigma - 0.5)*sum(Elog_phi) - sum(Ephi.*(0.5*diag(Ealpha_alphaT)+hyper.varsigma)) - ...
       (N+1)*log(gamma(hyper.sigma+0.5)));
%Varphi
bound = bound -K*log(gamma(hyper.tau)) + K*hyper.tau*log(hyper.nu) + ...
       (1-hyper.tau)*sum(Elog_varphi) - hyper.nu * sum(Evarphi);
bound = bound - ( (hyper.tau - 1)*sum(Elog_varphi) + ...
       sum((Evarphi-1).*Elog_beta) - hyper.nu*sum(Evarphi) - log(Zphi) );    
%Beta
bound = bound + Elog_gamma_sum_varphi - sum(Elog_gamma_varphi) + sum((Evarphi-1).*Elog_beta);
bound = bound - ( -(1/2)*sum(sum(Ebeta_betaT.*Omega)) + ...
       0.5*Ebeta'*b + sum((Evarphi - 1).*Elog_beta) - log(Zbeta));
       
%Likelihood
bound = bound + sum(log(1./(1+exp(-Xi)))) + sum(Xi.^2) - sum(Xi./2) + 0.5 * sum(t'.*(Ebeta'*Z)) - ...
    0.5*sum(sum(Omega.*Ebeta_betaT));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ebeta, Elog_beta, Ebeta_betaT, SumW, Offset, exflag] = sample_beta(...
           NoS, Ebeta, Elog_beta, Ebeta_betaT, SumW, Evarphi, Omega, b, Offset)
    %Firstly draw samples from Dirichlet with pars Evarphi
    K = length(Ebeta);
    G = repmat(0,K,NoS);
    for k = 1:K
        %G(k,:) = gamrnd(Evarphi(k),1,1,NoS);
        if Evarphi(k)>1
          G(k,:) = randraw('gamma',[0,1,Evarphi(k)],NoS)';
        else
          G(k,:) = randraw('gamma',[0,1,Evarphi(k)+1],NoS)';
          G(k,:) = G(k,:).*(rand(1,NoS).^(1/Evarphi(k)));
        end
    end
    G = G./repmat(sum(G,1),K,1);
    prob = find(sum(G,1)==0);
    G(:,prob) = [];
    NoS = NoS - length(prob);
    G(find(G<1e-10)) = 1e-10;
    Q = -(1/2)*(diag(G'*Omega*G)-G'*b)';
    if Offset == 0
        Offset = max(Q);
    end
    W = exp(Q - Offset);
    if any(isinf(W))
        I = find(isinf(W));
        exflag = -1;
        return
    end
    Ebeta = (sum(G.*repmat(W,K,1),2)+Ebeta*SumW)/(sum(W) + SumW);
    Elog_beta = (sum(log(G).*repmat(W,K,1),2)+Elog_beta*SumW)/(sum(W) + SumW);
    Ebeta_betaT = (((G.*repmat(W,K,1))*G') + Ebeta_betaT*SumW)/(sum(W) + SumW);
    exflag = 1;
    
    SumW = SumW + sum(W);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Evarphi, Elog_varphi, Elog_gamma_varphi, Elog_gamma_sum_varphi, SumW, Offset] = sample_varphi(...
           NoS, Evarphi, Elog_varphi, Elog_gamma_varphi, Elog_gamma_sum_varphi,...
           SumW,tau,nu,Elog_beta,Offset)
    %Firstly, take NoS samples from Gamma(tau,nu)
    %NOTE matlab defines Gamma(tau,1/nu)
    K = length(Evarphi);
    %B = gamrnd(tau,1/nu,K,NoS);
    for k = 1:K
      if tau > 1
        B(k,:) = randraw('gamma',[0,1/nu,tau],NoS)';
      else
        B(k,:) = randraw('gamma',[0,1/nu,1+tau],NoS)';
        B(k,:) = B(k,:).*(rand(1,NoS).^(1/tau));
      end
    end
    %Calculate W(varphi)
    G = gamma(sum(B,1))./prod(gamma(B),1);
    T = sum((B-1).*repmat(Elog_beta,1,NoS),1);
    %if Offset == 0
    %    Offset = max(T);
    %    NewOffset = max(T);
    %else
    %    NewOffset = max(T);
    %end
    NewOffset = 0;
    W = G.*exp(T-max(T));
    chg = exp(NewOffset - Offset);
    Evarphi = (sum(B.*repmat(W,K,1),2)*chg+Evarphi*SumW)./(chg*sum(W)+SumW);
    Elog_varphi = (sum(safe_log(B).*repmat(W,K,1),2)*chg+Elog_varphi*SumW)./(chg*sum(W)+SumW);
    Elog_gamma_varphi = (sum(log(gamma(B)).*repmat(W,K,1),2)*chg+Elog_gamma_varphi*SumW)./(chg*sum(W)+SumW);
    Elog_gamma_sum_varphi = (sum(log(gamma(sum(B,1))).*W,2)*chg+Elog_gamma_sum_varphi*SumW)./(chg*sum(W)+SumW);
    SumW = SumW + chg*sum(W);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = Lam(X)

o = (1./(4*X)).*tanh(0.5*X);
o(find(o>1e10)) = 1e10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = safe_log(X)

X(find(X<1e-10)) = 1e-100;
o = log(X);