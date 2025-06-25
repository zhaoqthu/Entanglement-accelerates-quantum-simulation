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
NC= C*B- B*C; 

d10= opnorm(Array(C),2);
d11= opnorm(Array(C1),2);


ComAB=2*(hx*hy*n+hy*2*(n-1));
ComBBA=4*hx*hy^2*n+8*hy^2*2*(n-1);
ComAAB=4*hx^2*hy*n+8*hy*(n-1)+8*hy*hx*2*(n-1)+8*hy*(n-2);
normA=hx*n+n-1;
normB=hy*n;

ave1=sqrt(4*(n-1)*hy^2*2+4*hy^2*hx^2*n)*dt^2/2;
wor1=2*(2*hy*(n-1)+hy*hx*n)*dt^2/2;

ave21=sqrt(16*hx^2*hy^4*n+16*4*hy^4*2*(n-1))*dt^3/12;
ave22=sqrt(16*(hx^2+2)^2*hy^2*(n)+16*4*2*hx^2*hy^2*(n-1)+16*4*hy^2*(n-2))*dt^3/24;

ave2=ave21+ave22
println(ave2)


worst2(r)=r*((4*hx*hy^2*n+8*hy^2*2*(n-1))*(t/r)^3/12 + (4*hx^2*hy*(n)+4*hy*n*2+8*hx*hy*(n-1)*2+8*hy*(n-2))*(t/r)^3/24);



cp=4;
#number of check point
theo_r=zeros(cp)
for j in 1:cp
  println(j)
  theo_r[j]=total_step_entropycheck2(A, B, n, hx, hy, t, j, 0.00001);
end
#according to different estimation method can choose entropycheck1 or entropycheck2
println(theo_r)

Stepw=zeros(1);
Stepw= total_theo_step_wor2(n, hx, hy, t, 0.00001)
println("WOR:", Stepw)
Stepa=zeros(1);
Stepa= total_theo_step_ave2(n, hx, hy, t, 0.00001)
println("ave:", Stepa)
#end





