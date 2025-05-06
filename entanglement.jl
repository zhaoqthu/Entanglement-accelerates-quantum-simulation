
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
ma=100;
rm=1;   #different magnetic
n=12;

#tem=n;
d00=zeros(1);
d01=zeros(1);
da0=zeros(1);
da1=zeros(1);
da20=zeros(1);
da21=zeros(1);

d10=zeros(ma);
d11=zeros(ma);
d20=zeros(ma);
d21=zeros(ma);


v1=zeros(ma);
v2=zeros(ma);
v3=zeros(ma);
v4=zeros(ma);
u1=zeros(ma);
u2=zeros(ma);
u3=zeros(ma);
u4=zeros(ma);
v21=zeros(ma);
d30=zeros(1);
d31=zeros(1);


#for j1 in 1:rm
t=10;
dt=t/ma;
dim = 2^n;
d = dim;
d1= sqrt(d);
#se=HaarKet{2}(d);
psi1=ZeroInp(n);
#a0=spzeros(ComplexF64, 1, d);
#b0=spzeros(ComplexF64, 1, d);
A = spzeros(ComplexF64, d, d);
B = spzeros(ComplexF64, d, d);
A1 = spzeros(ComplexF64, d, d);
B1 = spzeros(ComplexF64, d, d);

#tol=1e-3;.

init_QIMF!(A, B, n, 0.8090, 0.9045, 1.0);
init_QIMF!(A1, B1, n, 0.1, 0.9045, 1.0);
H=A+B;
H1=A1+B1;
#C= A * B - B * A;
#C1=A1 * B1 - B1 * A1;
C=exp(-1im .* dt*Array(H))-exp(-1im .* dt*Array(A))*exp(-1im .* dt*Array(B))
C1=exp(-1im .* dt*Array(H1))-exp(-1im .* dt*Array(A1))*exp(-1im .* dt*Array(B1))

C20=exp(-1im .* dt*Array(H))-exp(-1im .* dt/2*Array(A))*exp(-1im .* dt*Array(B))*exp(-1im .* dt/2*Array(A))
C21=exp(-1im .* dt*Array(H1))-exp(-1im .* dt/2*Array(A1))*exp(-1im .* dt*Array(B1))*exp(-1im .* dt/2*Array(A1))


#com_ABB= C*B- B*C;
#com_AAB= C*A- A*C;
##com_ABB1= C1*B1- B1*C1;
#com_AAB1= C1*A1- A1*C1;



#d00= opnorm(Array(C))*dt^2/2;
d00= opnorm(C);
println(d00)
d20= opnorm(C20);
println(d20)

#d01= opnorm(Array(C1))*dt^2/2;
d01= opnorm(C1);
println(d01)
d21= opnorm(C21);
println(d21)
#dw20= dt^3 * opnorm(Array(com_AAB))/24+ dt^3 * opnorm(Array(com_ABB))/12;
#println(dw20)
#dw21= dt^3 * opnorm(Array(com_AAB1))/24+ dt^3 * opnorm(Array(com_ABB1))/12;
#println(dw21)



#da0 = sqrt(sp_trace(C * C')/d)*dt^2/2;
#da1=  sqrt(sp_trace(C1 * C1')/d)*dt^2/2;


da0 = sqrt(tr(C * C')/d);
da1=  sqrt(tr(C1 * C1')/d);
println(da0)
println(da1)

da20 = sqrt(tr(C20 * C20')/d);
da21=  sqrt(tr(C21 * C21')/d);
println(da20)
println(da21)

#da20 = sqrt(sp_trace(com_ABB* com_ABB')/d)*dt^3/12+ sqrt(sp_trace(com_AAB * com_AAB')/d)*dt^3/24;
#da21=  sqrt(sp_trace(com_ABB1* com_ABB1')/d)*dt^3/12+ sqrt(sp_trace(com_AAB1 * com_AAB1')/d)*dt^3/24;
#println(da20)
#rintln(da21)

#f1=Matrix(C*C');
#d30=sqrt(tr(f1)/d)
#d31=sqrt(trace_sparse_complex_matrix(C*C')/d)


#dt=0.1;
for j in 1:ma

println(j)
a0=expmv(-1im .* t*(j)/100, H, psi1);
b0=expmv(-1im .* t*(j)/100, H1, psi1);

#partial_trace(a0*a0',n,[2])
a10=expmv(-1im .* dt, H, a0);
a1a=expmv(-1im .* dt, A, a0);
a11=expmv(-1im .* dt, B, a1a);
d10[j]= norm(a10-a11);
#println(d10)
a20=expmv(-1im .* dt/2, A, a0);
a21=expmv(-1im .* dt, B, a20);
a22=expmv(-1im .* dt/2, A, a21);
d20[j]= norm(a10-a22);

b10=expmv(-1im .* dt, H1, b0);
b1a=expmv(-1im .* dt, A1, b0);
b11=expmv(-1im .* dt, B1, b1a);
d11[j]= norm(b10-b11);

b20=expmv(-1im .* dt/2, A1, b0);
b21=expmv(-1im .* dt, B1, b20);
b22=expmv(-1im .* dt/2, A1, b21);
d21[j]= norm(b10-b22);

rhoa4=TrX4R(a0*a0',1,2,3,4,n);

rhoa3=TrX4R(a0*a0',1,2,3,3,n);
rhoa2=TrX4R(a0*a0',1,2,2,2,n);
rhoa1=TrX4R(a0*a0',n,n,n,n,n);

sig4=TrX4R(b0*b0',1,2,3,4,n);

sig3=TrX4R(b0*b0',1,2,3,3,n);
sig2=TrX4R(b0*b0',1,2,2,2,n);
sig1=TrX4R(b0*b0',n,n,n,n,n);

#println(size(rhoa4))
#rho1=TrX4R(b0*b0',1,2,3,4,n);
v1[j]=von_neumann_entropy(rhoa1);
#println(rhoa1)
#println(size(rhoa2))
v2[j]=von_neumann_entropy(rhoa2);
v3[j]=von_neumann_entropy(rhoa3);
v4[j]=von_neumann_entropy(rhoa4);

u1[j]=von_neumann_entropy(sig1);
#println(v4)
#println(size(rhoa2))
u2[j]=von_neumann_entropy(sig2);
u3[j]=von_neumann_entropy(sig3);
u4[j]=von_neumann_entropy(sig4);
end



x=10/ma:10/ma:10;
plot!(x,d10)
plot!(x,d11)
plot!(x,d20)
plot!(x,d21)
plot!(x,v3)
plot!(x,v4)
plot!(x,v2)
plot!(x,v1)


plot!(x,u3)
plot!(x,u4)
plot!(x,u2)
plot!(x,u1)
#end
#println(B)
#init_powerlaw!(A,B,C,n,0)
#RandInp(n, 0)
#RandInp(n, 1)
#  #println(C)
#se=HaarKet{2}(d);
#psi1=rand(se);
#psi1=ket(3,d)
#H=A+B;
#

#ep2[j]=trotteran6(A,B,t,,40);

#mean(r1[j,:])


#@time epp[j]=trotteran(A,B,t,r,20);

#trotteran4
#epi[j]= triangle_error(A,B,t,r);
#println(ep[j])
#7e+4,9e+4
#@time r[j]=search_empir_numR(A,B,tem,tol,100,200);
#r2[j]=search_trotter_num_triangle3(A,B,C,tem,tol);
#println(r[j])
#println(epp[j])
#end


#H=A+B;
#r1=search_trotter_num_triangle(A, B, t, tol);
#r2=search_trotter_num_inteference(A, B, t, tol);

#@time   =trottervec(A,B,t,)
#@time a4=trotteran(A, B, t, 200, 20)
#@time a5=trotteran(A, B, t, 27400, 20)
#println(a5)
#@time r3=search_empir_numR(A,B,t,tol,1700,2400);
#@time r3=search_empir_numR(A,B,t,tol,100,0.5e+4);
