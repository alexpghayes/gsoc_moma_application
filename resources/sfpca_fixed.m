%%%%%%%%%%%%%%%%%%
%function to compute the rank-one SFPCA solution
%L1 sparse penalties
%for fixed regularization params
%%%%%%%%%%%%%%%%%%%%%%%%
%inputs:
%x: n x p data matrix 
%lamu: sparsity parameter for u
%lamv: sparsity parameter for v
%alphau: smoothness parameter for u
%alphav: smoothness parameter for v
%Omegu: n x n positive semi-definite matrix for roughness penalty
%on u
%Omegv: p x p positive semi-definite matrix for roughness penalty
%on v
%startu: n  vector of starting values for U;  if startu=0, then
%algorithm initialized to the rank-one SVD solution
%startv: p vector of starting values for V;  if startv=0, then
%algorithm initialized to the rank-one SVD solution
%posu: non-negativity indicator - posu = 1 imposes non-negative
%constraints on u, posu = 0 otherwise
%posv: non-negativity indicator - posv = 1 imposes non-negative
%constraints on v, posv = 0 otherwise
%maxit: maximum number of alternating regressions steps
%%%%%%%%%%%%%%%%
%outputs:
%U: n x 1 left SFPC
%V p x 1 right SFPC
%d: the associated singular value
%Xhat: the deflated residual matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[U,V,d,Xhat] = sfpca_fixed(x,lamu,lamv,alphau,alphav,Omegu,Omegv,startu,startv,posu,posv,maxit)
[n,p] = size(x);
Su = speye(n) + n*alphau*Omegu;
Sv = speye(p) + p*alphav*Omegv;
Lu = max(eig(Su)) + .01;
Lv = max(eig(Sv)) + .01;
thr = 1e-6;  
Xhat = x; 
if sum(startu)==0
    [u,d,v] = svds(x,1);
else
    u = startu;
    v = startv;
end
indo = 1; iter = 0;
while indo>thr & iter<maxit
    oldu = u; oldv = v;
    indu = 1;
    while indu>thr
        oldui = u;
        utild = u + (Xhat*v - Su*u)/Lu;
        u = soft_thr(utild,lamu/Lu,posu);
        if norm(u)>0
            u  = u/sqrt(u'*Su*u);
        else
            u = zeros(n,1);
        end            
        indu = norm(u - oldui)/norm(oldui);
    end
    indv = 1;
    while indv>thr
        oldvi = v;
        vtild = v + (Xhat'*u - Sv*v)/Lv;
        v = soft_thr(vtild,lamv/Lv,posv);
        if norm(v)>0
            v  = v/sqrt(v'*Sv*v);
        else
            v = zeros(p,1);
        end            
        indv = norm(v - oldvi)/norm(oldvi);
    end
    indo = norm(oldu - u)/norm(oldu) + norm(oldv - v)/norm(oldv);
    iter = iter + 1
end
U = u/norm(u);
V = v/norm(v);
d = U'*Xhat*V;
Xhat = Xhat - d*U*V';

load fisheriris.mat
sfpca_fixed(x = meas, lamu = 0,lamv = 0,alphau = 0,alphav = 0, %
            Omegu = ,Omegv,startu,startv,posu,posv,maxit)
