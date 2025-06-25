using LinearAlgebra
using SparseArrays
using ExpmV
using Plots
using QuantumInformation
using Statistics
#using Yao, BitBasis
#using LazySets
#include("Qudos.jl")
include("simulation_entanglement.jl")

using Main.simulation_entanglement
ma=50;
rm=1;   #different magnetic
n=12;
dt=0.001;

hx=0.8090;
#hx=0.0;
hy=0.9045;
J=1.0
#tem=n;
d10=zeros(1);
d11=zeros(1);
d20=zeros(ma);
d21=zeros(ma);
v20=zeros(ma);
e1DS=zeros(ma);
e2DS1=zeros(ma);
e2DS2=zeros(ma);
v21=zeros(ma);
d30=zeros(1);
d31=zeros(1);
da1=zeros(ma);
dw1=zeros(ma);
TEnB1=zeros(ma);

#for j1 in 1:rm
t=n;
dim = 2^n;
d = dim;
d1= sqrt(d);
#se=HaarKet{2}(d);
psi1=ZeroInp(n);

A = spzeros(ComplexF64, d, d);
B = spzeros(ComplexF64, d, d);
A1 = spzeros(ComplexF64, d, d);
B1 = spzeros(ComplexF64, d, d);

#tol=1e-3;.

init_QIMF!(A, B, n, hx, hy, J);
#init_QIMF!(A1, B1, n, 0.1, 0.9045, 1.0);
H=A+B;
H1=A1+B1;
C= A * B - B * A;
C1=A1 * B1 - B1 * A1;
NC1= C*B- B*C; 
NC2= C*A- A*C;
NC11= NC1*B- B*NC1;
NC12=NC1*A- A*NC1;
NC21= NC2*B- B*NC2;
NC22= NC2*A- A*NC2;
d10= opnorm(Array(C),2);
d11= opnorm(Array(C1),2);

#ComAB=d10;
#normA=opnorm(Array(A),2);
#normB=opnorm(Array(B),2);

#ComBBA=opnorm(Array(NC1),2);

#ComAAB=opnorm(Array(NC2),2);

#ComBBBA=opnorm(Array(NC11),2);

#ComABBA=opnorm(Array(NC12),2);
#ComBAAB=opnorm(Array(NC21),2);
#ComAAAB=opnorm(Array(NC22),2);


ComAB=2*(hx*hy*n+hy*2*(n-1));
ComBBA=4*hx*hy^2*n+8*hy^2*2*(n-1);
ComAAB=4*hx^2*hy*n+8*hy*(n-1)+8*hy*hx*2*(n-1)+8*hy*(n-2);



normA=hx*n+n-1;
normB=hy*n;

ComBBBA=ComBBA*2*normB;

ComABBA=ComBBA*2*normA;
ComBAAB=ComAAB*2*normB;
ComAAAB=ComAAB*2*normA;

ave1=sqrt(4*(n-1)*hy^2*2+4*hy^2*hx^2*n)*dt^2/2;
wor1=2*(2*hy*(n-1)+hy*hx*n)*dt^2/2;

ave21=sqrt(16*hx^2*hy^4*n+16*4*hy^4*2*(n-1))*dt^3/12;
ave22=sqrt(16*(hx^2+2)^2*hy^2*(n)+16*4*2*hx^2*hy^2*(n-1)+16*4*hy^2*(n-2))*dt^3/24;

#ave21=sqrt(sp_trace( NC1* NC1')/(2^n)) * dt^3/12;
#ave22=sqrt(sp_trace( NC2* NC2')/(2^n)) * dt^3/24;

ave2=ave21+ave22
println(ave2)
wor2=(4*hx*hy^2*n+8*hy^2*2*(n-1))*dt^3/12 + (4*hx^2*hy*(n)+4*hy*n*2+8*hx*hy*(n-1)*2+8*hy*(n-2))*dt^3/24;
#wor2=ComBBA*dt^3/12+ ComAAB*dt^3/24;
println(wor2)
#wor1=d10*dt^2/2;
#f1=Matrix(C*C');
#d30=sqrt(tr(f1)/d)
#d31=sqrt(trace_sparse_complex_matrix(C*C')/d)


#dt=0.1;
for j in 1:50

a0=expmv(-1im .*(dt*j*200), H, psi1);
#b0=expmv(-1im .* t*j/100, H1, psi1);
#println(a0)
#partial_trace(a0*a0',n,[2])
#a10=expmv(-1im .* dt, H, a0);
#a1a=expmv(-1im .* dt, A, a0);
#a11=expmv(-1im .* dt, B, a1a);

#b10=expmv(-1im .* dt, H1, b0);
#b1a=expmv(-1im .* dt, A1, b0);
#b11=expmv(-1im .* dt, B1, b1a);
#d20[j]= norm(a10-a11);
#v20println(d20)
#d21[j]= norm(b10-b11);
#d20[j]= norm(C*a0);
#d21[j]= norm(C1*b0);
#println(size(a0*a0'))
#rho0=TrX4R(a0*a0',1,1,1,1,n);
#rho=a0*a0';
da1[j]=ave1;
dw1[j]=wor1;

#println(tr(rho))
#entS = 0
#ent =0
#for a1 in 1:(n-1)
#  a2 = a1 + 1
#  for a3 in 1:(n-1)
#      a4 = a3 + 1
#      b = sort([a1, a2, a3, a4])
    #println(b)
 #     rho4 = TrX4R(rho, b[1], b[2], b[3], b[4], n)
      #rho5=  TrX4R(rho, 1, 1, 2, 2, n)
     # println(tr(rho5))
  #    d1, _ = size(rho4)
   #   Dis = rho4 - I(d1) / d1
      #ent= norm(Dis)
    #  ent = tr(sqrt(Dis * Dis'))
    #  entS = entS + ent
      #println(ent)
      #println(entS)
  #end
#end
#v20[j]=entS;
#println(v20)
#println(entS*4*hy^2*4)
#e1DS[j]=  entS_DisSD2(a0*a0', n,4*hy^2,4 ) + entS_DisS21(a0*a0', n,4*hy^2*hx,4 )+ entS_DisS11(a0*a0', n,4*hy^2*hx^2,1 )



D22=entS_DisSD2(a0*a0', n,16*4*hy^4,2*2 );
D21=entS_DisS21(a0*a0', n,64*hx*hy^4,2*2 );
D11=entS_DisS11(a0*a0', n,16*hy^4*hx^2,1 );
#println(D22)
#println(D21)
#println(D11)


e2DS1[j]=D22+D21+D11;
#println(e2DS1[j])
#println(a0)
D33=entS_DisSD3(a0*a0', n, (8*hy)^2,1 );
D32=entS_DisS32(a0*a0', n, 64*hx*hy^2,2*2 );
D31=entS_DisS31(a0*a0', n, 32*(hx^2+2)*hy^2,2 );
D33=real(D33);
D32=real(D32);
D31=real(D31);
#println(D33)
#println(D32)
#println(D31)

e2DS2[j]=D33+D32 +D31 + D22/hy^2*hx^2  +  D21/(hx*hy^4)*(hx^2*hy+2*hy)*hx*hy+ D11/(hy^4*hx^2)*(hx^2*hy+2*hy)^2 ;  



println(e2DS2[j])

#TEnB1[j]=sqrt(ave21^2+ dt^6/144*e2DS1[j]+(5*dt^7/448*(2*normB+normA)+dt^8/120*(normB^2+normA^2/4))*ComBBA^2)+ sqrt( ave22^2+ dt^6/576*e2DS2[j]+(5*dt^7/1792*(2*normB+normA)+ dt^8/480*(normB^2+normA^2/4))*ComAAB^2 )+ dt^4/48*ComBBA+ dt^4/192*ComAAB;
TEnB1[j]=sqrt(ave21^2+ dt^6/144*e2DS1[j])+ sqrt( ave22^2+ dt^6/576*e2DS2[j])+ dt^4*ComABBA/32+ dt^4*ComBBBA/12+ dt^4*ComBAAB/32 + dt^4*ComAAAB/48 ;

#TEnB1[j]=sqrt(ave1^2+ dt^4/4*e1DS[j]+ dt^5/5*3*normA*ComAB^2+ dt^6/72*ComAB^3)+ dt^3/6*ComBBA;

#TEnB1[j]=sqrt(ave1^2+ dt^4/4*entS*4*hy^2*4+ dt^5/5*3*normA*ComAB^2+ dt^6/72*ComAB^3)+ dt^3/6*ComBBA;
println(TEnB1)
#M = zeros(Complex{Float64},  2^(d-a2+a1-1), 2^(d-a2+a1-1))
#for i in 1:2^(a2-a1 + 1)
 #   mi=zeros(Complex{Float64},1, 2^(a2-a1 + 1))
  #  left=I(2^(a1-1));
  # mi[i]=1
  # right=I(2^(d-a2));
  # V=kron(left,mi,right);
  # M +=V*rho*V'
  # println(M)
  # printlnvonneumann_entropy(M)
#end
#println(M)
#rho1=ptrace(a1*a1', [2^5, 2^5], [2]);


#v21[j]=vonneumann_entropy(rho1);

end



x=0.1:0.1:10;
plot!(x,v20)
plot!(x,TEnB1)
plot!(x,dw1)
#plot!(x,d21)

#end

#v20*4*hy^2*4
#sqrt(ave1^2+ dt^4/4*e1DS[j]+ dt^5/5*3*normA*ComAB^2+ dt^6/72*ComAB^3)+ dt^3/6*ComBBA;