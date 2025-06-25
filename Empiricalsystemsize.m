clc;
clear;

nmax = 10;
epsilon = 1e-3;
em_trotter_num = zeros(1,nmax-2);
th_trotter_num = zeros(1,nmax-2);

em_trotter_numIn = zeros(1,nmax-2);
em_trotter_numRan = zeros(1,nmax-2);


for n = 3:nmax
[A,B]=QIMF(n,0.8090, 0.9045, 1);

  fprintf('n = %d\n', n);
 
    init_guess = n * 1e10;
em_trotter_numIn(n-2)=empirical_trotter_numberIn(A, B, epsilon, init_guess, 1)
em_trotter_numRan(n-2)=empirical_trotter_numberRan(A, B, epsilon, init_guess, 20)
   em_trotter_num(n-2)=empirical_trotter_number2(A, B, epsilon, init_guess)
%th_trotter_num(n-2)=theoretial_trotter_numberIn(A, B, epsilon, 2*em_trotter_num(n-2), 1)
end    


nmax =10
loglog(3:nmax, em_trotter_numIn(1:nmax-2),'s','MarkerEdgeColor',[0.4940    0.1840    0.5560 ],'MarkerFaceColor',[0.4940    0.1840    0.5560] )
hold on
loglog(3:nmax, em_trotter_numRan(1:nmax-2),'s','MarkerEdgeColor',[0.4660    0.6740    0.1880],'MarkerFaceColor',[0.4660    0.6740    0.1880] )
hold on

loglog(3:nmax, em_trotter_num(1:nmax-2),'s', 'MarkerEdgeColor',[100,149,237]/256,'MarkerFaceColor',  [100,149,237]/256)
hold on


p_In = polyfit(log(3:nmax), log(em_trotter_numIn(1:nmax-2)),1) %random impirical
p_Ran = polyfit(log(3:nmax), log(em_trotter_numRan(1:nmax-2)),1)
p_o= polyfit(log(3:nmax), log(em_trotter_num(1:nmax-2)),1)
ncoord = 3:12
loglog(ncoord, ncoord.^p_In(1) .*exp(p_In(2)),'k--');
loglog(ncoord, ncoord.^p_Ran(1) .*exp(p_Ran(2)),'k--');
loglog(ncoord, ncoord.^p_o(1) .*exp(p_o(2)),'k--');
%10^p_em(1)*exp(p_em(2))
text(12,12^p_In(1) .*exp(p_In(2)*1.01),'$\mathcal O(n^{  1.4426  })$','interpreter','latex','FontSize',13);
text(12,12^p_Ran(1) .*exp(p_Ran(2)/1.01),'$\mathcal O(n^{ 1.4354})$','interpreter','latex','FontSize',13);
text(12,12^p_o(1) .*exp(p_o(2)),'$\mathcal O(n^{ 1.7797})$','interpreter','latex','FontSize',13);

xlabel("Number of qubits $n$",'interpreter','latex','fontsize',18)
ylabel("$r$",'interpreter','latex','fontsize',18)
%title('Power law $\alpha=4$ PF1','interpreter','latex')
title('QIMF $h_x=0.8$ PF1','interpreter','latex','fontsize',18)
% %loglog(ncoord, ncoord.^p_th(1) .*exp(p_th(2)),'b--');
legend('Location','northwest')
legend('Empirical with input','Empirical with random inputs','Empirical (spectral norm)','Extrapolation curve')


%% Compute empirical trotter number
%  nmax = 5;
% num_trials = 1;
% num_samples = 20;
% epsilon = 1e-3;
% em_trotter_num = zeros(1,nmax-2);
% 
% for n = 3:nmax
%     fprintf('n = %d\n', n);
%     holder = zeros(1,num_trials);
%     init_guess = n * 1e3;
%     
%     for counter = 1:num_trials
%         [A,B] = heisenberg(n);
%         holder(counter) = empirical_trotter_number(A, B, epsilon, init_guess, num_samples);
%     end
%     
%     em_trotter_num(n-2) = mean(holder);
%     fprintf("\n")
% end
% % 
% %% save data to MAT file
% %save('em_trotter_number.mat','em_trotter_num')
% %% Compute theoretical trotter number
% nmax = 12;
% epsilon = 1e-3;
% num_trials = 5;
% th_trotter_num = zeros(1,nmax-2);
% 
% for n = 3:nmax
%     fprintf('n = %d\n', n);
%     init_guess = n * 1e4;
%     holder2 = zeros(1,num_trials);
%     
%     for counter = 1:num_trials
%         [A,B] = heisenberg(n);
%         holder2(counter) = theoretical_trotter_number(A, B, epsilon, init_guess);
%     end
%     th_trotter_num(n-2) = mean(holder2);
%     fprintf("\n")
% end
% 
% %% linear regression
% nmax=7;
% p_em = polyfit(log(3:nmax), log(em_trotter_num(1:nmax-2)),1);
% p_th = polyfit(log(3:12), log(th_trotter_num),1);
% ncoord = 3:100;
% 
% %% make plots
% 
% figure;
% loglog(3:11, em_trotter_num(1:9),'r-s');
% hold on
% loglog(ncoord, ncoord.^p_em(1) .*exp(p_em(2)),'r--');
% loglog(3:12,th_trotter_num,'b-o');
% loglog(ncoord, ncoord.^p_th(1) .*exp(p_th(2)),'b--');
% legend('empirical','extrapolation-empitical','our bounds','extrapolation-bound')

%%
function [A,B] = heisenberg(n,k)
X = sparse([0 1;1 0]);
Y = sparse([0 -1i;1i 0]);
Z = sparse([1 0;0 -1]);
XX = kron(X,X);
YY = real(kron(Y,Y));
ZZ = kron(Z,Z);

A = sparse(2^n,2^n);
B = sparse(2^n,2^n);

for j = 1:floor(n/2)
    A = A + kron(speye(2^(2*j-2)), kron(XX, speye(2^(n - 2*j)) ) );
    A = A + kron(speye(2^(2*j-2)), kron(YY, speye(2^(n - 2*j)) ) );
    A = A + kron(speye(2^(2*j-2)), kron(ZZ, speye(2^(n - 2*j)) ) );
    A = A + k*(2*rand(1)-1) .* kron(speye(2^(2*j-2)), kron(Z, speye(2^(n - 2*j + 1)) ) );
end

for j = 1:(ceil(n/2)-1)
    B = B + kron(speye(2^(2*j-1)), kron(XX, speye(2^(n - 2*j-1)) ) );
    B = B + kron(speye(2^(2*j-1)), kron(YY, speye(2^(n - 2*j-1)) ) );
    B = B + kron(speye(2^(2*j-1)), kron(ZZ, speye(2^(n - 2*j-1)) ) );
    B = B + k*(2*rand(1)-1) .* kron(speye(2^(2*j-1)), kron(Z, speye(2^(n - 2*j)) ) );
end
end


function [A,B] = QIMF(n,hx, hy,J)
X = sparse([0 1;1 0]);
Y = sparse([0 -1i;1i 0]);
Z = sparse([1 0;0 -1]);
XX = kron(X,X);
%YY = real(kron(Y,Y));
%ZZ = kron(Z,Z);

A = sparse(2^n,2^n);
B = sparse(2^n,2^n);

for j = 1: n-1
    A = A + hx.*kron(speye(2^(j-1)), kron(X, speye(2^(n - j)) ) );
    A = A + J.* kron(speye(2^(j-1)), kron(XX, speye(2^(n -1- j)) ) );
end
A= A+ hx.*kron(speye(2^(n-1)), X);

for j = 1:n 
    B = B + hy.*kron(speye(2^(j-1)), kron(Y, speye(2^(n - j)) ) );
end
end



function [A,B,C] = powerlaw(n,alpha)
X = sparse([0 1;1 0]);
Y = sparse([0 -1i;1i 0]);
Z = sparse([1 0;0 -1]);
%XX = kron(X,X);
%YY = real(kron(Y,Y));
%ZZ = kron(Z,Z);

A = sparse(2^n,2^n);
B = sparse(2^n,2^n);
C = sparse(2^n,2^n);

for j = 1:n-1
    for i= j+1:n
      XX= kron(X, kron(speye(2^(i-j-1)),X))./(i-j)^alpha; 
      YY= kron(Y, kron(speye(2^(i-j-1)),Y))./(i-j)^alpha; 
      ZZ= kron(Z, kron(speye(2^(i-j-1)),Z))./(i-j)^alpha; 
    A = A + kron(speye(2^(j-1)), kron(XX, speye(2^(n - i))));
    B = B + kron(speye(2^(j-1)), kron(YY, speye(2^(n - i))));
    C = C + kron(speye(2^(j-1)), kron(ZZ, speye(2^(n - i))));
    end
end

end



function err = theoretical_error2(A, B, r, psi00)
%%
[d,~] = size(A);
n = log2(d);
H=A+B;
t=n;
hx=0.8090;
hy=0.9045;
dt=t/r;
U2 = expm(-1i*dt*A/2) * expm(-1i*dt*B)* expm(-1i*dt*A/2);
U0 = expm(-1i*n*H);
Com=A*B-B*A;
ComBBAP=B*Com-Com*B;
ComAABP=A*Com-Com*A;
E2=  U2-U0;
wor=norm(full(E2));

ave21=norm(ComBBAP,'fro')./sqrt(2^n)*dt^3/12;

ave22=norm(ComAABP,'fro')./sqrt(2^n)*dt^3/24;

ComAB=norm(full(Com),2);
ComBBA=norm(full(ComBBAP),2);
ComAAB=norm(full(ComAABP),2);
normA=norm(full(A),2);
normB=norm(full(B),2);
tterror=0;

for j=1:1:r
    j
psi=U0^j*psi00;
D22=DisSD2(psi*psi', n, 16*4*hy^4,2*2);
D21= DisS21(psi*psi', n,  64*hx*hy^4,2*2);
D11=DisS11(psi*psi', n, 16*hy^4*hx^2,1);
e2DS1=  D22+  + D21+  D11;

D33=DisSD3(psi*psi', n, (8*hy)^2,1 );
D32=DisS32(psi*psi', n, 64*hx*hy^2,2*2 );
D31=DisS31(psi*psi', n, 32*(hx^2+2)*hy^2,2 );
e2DS2=  D33+D32 +D31+ D22/hy^2*hx^2 +  D21/(hx*hy^4)*(hx^2*hy+2*hy)*hx*hy+ D11/(hy^4*hx^2)*(hx^2*hy+2*hy)^2 ;   %+ DisS5(psi*psi', n, (4*hy)^2, 2+2*2)+ DisS32(psi*psi', n, (4*hy)^2, 2+2*2) +      DisSD2(psi*psi', n)+DisS3(psi*psi', n)+DisS2(psi*psi', n);

error=sqrt(ave21^2+ dt^6/144*e2DS1+(dt^8/120*(normB^2+normA^2/4)+dt^7*5/448*(2*normB+normA))*ComBBA^2)+ sqrt( ave22^2+ dt^6/576*e2DS2+ (dt^7*5/1792*(2*normB+normA)+dt^8/480*(normB^2+normA^2/4))*ComAAB^2 )+ dt^4/48*ComBBA+ dt^4/192*ComAAB;

error=min(error,wor);
tterror=tterror+error;

end


err=tterror

end






function err = theoretical_error(A, B, r)
%%
[d,~] = size(A);
n = log2(d);
%r = 1e6;
dt = n/r;
C = B*A - A*B;
step_err = 0.5 * dt^2 * sqrt(trace(C * C')/d);
err = r * step_err;
end

function err = empirical_error(A,B, r, psi)
%%
[d,~] = size(A);
n = log2(d);
dt = 1.5*n/r;

U0 = expm(-1i*dt*A) * expm(-1i*dt*B);
U = expm(-1i*n*(A+B));
err = norm(U0^r * psi - U * psi, 2);
end

function err = empirical_error2(A,B, r, psi)
%%
[d,~] = size(A);
n = log2(d);
dt = n/r;

U0 = expm(-1i*dt*A/2) * expm(-1i*dt*B)* expm(-1i*dt*A/2);
U = expm(-1i*n*(A+B));
err = norm(U0^r * psi - U * psi, 2);
end


function err = empirical_errorW(A,B, r)
%%
[d,~] = size(A);
n = log2(d);
dt = n/r;

U0 = expm(-1i*dt*A) * expm(-1i*dt*B);
U = expm(-1i*n*(A+B));
err = norm(U0^r - U , 2);
end

function err = empirical_errorW2(A,B, r)
%%
[d,~] = size(A);
n = log2(d);
t=n;
dt = t/r;

U0 = expm(-1i*dt*A/2) * expm(-1i*dt*B)*expm(-1i*dt*A/2) ;
U = expm(-1i*t*(A+B));
err = norm(U0^r - U , 2);
end




function rmin= empirical_errorWS(A,B, r)
%%
[d,~] = size(A);
n = log2(d);
dt = n/r;

U0 = expm(-1i*dt*A) * expm(-1i*dt*B);
U = expm(-1i*n*(A+B));
err = norm(U0^r - U , 2);
end

function rmin = binary_search(func, epsilon, init_r)
% func is a monotone decreasing function
left = 1;
right = init_r;

while func(right) > epsilon
    fprintf("Initial guess is too small! 10 times the initial guess!");
    right = 10 * init_r;
end

if func(left) <= epsilon
    rmin = 1;
else
    
    while left < right - 1
        mid = ceil((left + right)/2);
        if func(mid) > epsilon
            left = mid;
        else
            right = mid;
        end
     %   fprintf("left:%d, right:%d\n",left,right);
    end
    
    rmin = right;

end

end

function tn = empirical_trotter_numberIn(A, B, epsilon, init_r, num_samples)
[d,~] = size(A);
n = log2(d);
sample = zeros(1,num_samples);

for k = 1:num_samples
  %  psi0 = rand(d,1) + 1i * rand(d,1);
   % psi0 = psi0./norm(psi0);
   % psi0=binary_generate(d);
   ket0=[1,0]';
psi0=multi(n,ket0);
    func_em = @(r) empirical_error2(A, B, r, psi0);
    sample(k) = binary_search(func_em, epsilon, init_r);
end
tn = mean(sample);
end



function tn = theoretial_trotter_numberIn(A, B, epsilon, init_r, num_samples)
[d,~] = size(A);
n = log2(d);
sample = zeros(1,num_samples);

for k = 1:num_samples
  %  psi0 = rand(d,1) + 1i * rand(d,1);
   % psi0 = psi0./norm(psi0);
   % psi0=binary_generate(d);
   ket0=[1,0]';
psi0=multi(n,ket0);
    func_em = @(r) theoretical_error2(A, B, r, psi0);
    sample(k) = binary_search(func_em, epsilon, init_r);
end
tn = mean(sample);
end




function tn = empirical_trotter_numberRan(A, B, epsilon, init_r, num_samples)
[d,~] = size(A);
n = log2(d);
sample = zeros(1,num_samples);

for k = 1:num_samples
  %  psi0 = rand(d,1) + 1i * rand(d,1);
   % psi0 = psi0./norm(psi0);

    psi0=randominput(n);

    func_em = @(r) empirical_error2(A, B, r, psi0);
    sample(k) = binary_search(func_em, epsilon, init_r);
end
tn = mean(sample);
end

function [out] = binary_generate(nsize)

r = rand(nsize,1);
index = (r<0.5);
out = ones(size(r));
out(index) = zeros(length(r(index)),1);
end

function [out] = randominput(n)
ket0=[1,0]';
ket1=[0,1]';
out=1;
for i=1:1:n
r = rand(1,1);
if r<0.5
   out= kron(out,ket0);
else  
   out= kron(out,ket1);
end
end
end

function tn = empirical_trotter_number2(A, B, epsilon, init_r)
func_em = @(r) empirical_errorW2(A, B, r);
   tn= binary_search(func_em, epsilon, init_r);
end

function tn = empirical_trotter_number(A, B, epsilon, init_r)
func_em = @(r) empirical_errorW(A, B, r);
   tn= binary_search(func_em, epsilon, init_r);
end

function tn = theoretical_trotter_number(A, B, epsilon, init_r)
func_th = @(r) theoretical_error(A, B, r);
tn = binary_search(func_th, epsilon, init_r);
end
