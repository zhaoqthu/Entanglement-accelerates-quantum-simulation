module simulation_entanglement
export ZeroInp, entS_DisSD2, von_neumann_entropy, init_QIMF!,  RandInp, TrX4R, TrX6R, trace_sparse_complex_matrix, partial_trace

export entS_DisS21, entS_DisS11,sp_trace, search_trotter_num_triangle, triangle_errorSe,binary_search,search_trotter_num_triangleWor
export triangle_errorSeWor, commutator, search_empir_num, search_empir_numRan,trotter_errorwor2, trotter_errorave2, entS_DisSD3, entS_DisS32, entS_DisS31
export trotter_error_entropyPF2,  total_step_entropycheck, trotter_error_entropyPF2Para, total_theo_step_wor2, total_theo_step_ave2, total_step_entropycheck2
using LinearAlgebra
using SparseArrays
using ExpmV
using QuantumInformation
using Statistics
using Cubature
using Arpack





function trotter_error_entropyPF2(A, B,n, hx, hy, t, dt)  #trotter error in time t
    #d, ~ = size(A);
    H=A+B;
    #n = ceil(log2(d));
    
    psi1=ZeroInp(n);
  
   
    ComAB=2*(hx*hy*n+hy*2*(n-1));
    ComBBA=4*hx*hy^2*n+8*hy^2*2*(n-1);
    ComAAB=4*hx^2*hy*n+8*hy*(n-1)+8*hy*hx*2*(n-1)+8*hy*(n-2);
    normA=hx*n+n-1;
    normB=hy*n;
    ave21=sqrt(16*hx^2*hy^4*n+16*4*hy^4*2*(n-1))*dt^3/12;
    ave22=sqrt(16*(hx^2+2)^2*hy^2*(n)+16*4*2*hx^2*hy^2*(n-1)+16*4*hy^2*(n-2))*dt^3/24;
    a0=expmv(-1im .*t, H, psi1);
    D22=entS_DisSD2(a0*a0', n,16*4*hy^4,2*2 );
    D21=entS_DisS21(a0*a0', n,64*hx*hy^4,2*2 );
    D11=entS_DisS11(a0*a0', n,16*hy^4*hx^2,1 );
    e2DS1=D22+D21+D11;

    D33=entS_DisSD3(a0*a0', n, (8*hy)^2,1 );
    D32=entS_DisS32(a0*a0', n, 64*hx*hy^2,2*2 );
    D31=entS_DisS31(a0*a0', n, 32*(hx^2+2)*hy^2,2 );
    D33=real(D33);
    D32=real(D32);
    D31=real(D31);

    e2DS2=D33+D32 +D31 + D22/hy^2*hx^2  +  D21/(hx*hy^4)*(hx^2*hy+2*hy)*hx*hy+ D11/(hy^4*hx^2)*(hx^2*hy+2*hy)^2 ;  
    TEnB1=sqrt(ave21^2+ dt^6/144*e2DS1+(5*dt^7/448*(2*normB+normA)+dt^8/120*(normB^2+normA^2/4))*ComBBA^2)+ sqrt( ave22^2+ dt^6/576*e2DS2+(5*dt^7/1792*(2*normB+normA)+ dt^8/480*(normB^2+normA^2/4))*ComAAB^2 )+ dt^4/48*ComBBA+ dt^4/192*ComAAB;
    return TEnB1 
end


#calculate local RDM entropy information 
function trotter_error_entropyPF2Para(A, B,n, hx, hy, t)  #trotter error in time t
    #d, ~ = size(A);
    H=A+B;
    #n = ceil(log2(d));
    
    psi1=ZeroInp(n);
  
   
    a0=expmv(-1im .*t, H, psi1);
    rho=convert(Matrix{ComplexF64}, a0*a0');
    D22=entS_DisSD2(rho, n,16*4*hy^4,2*2 );
    D21=entS_DisS21(rho, n,64*hx*hy^4,2*2 );
    D11=entS_DisS11(rho, n,16*hy^4*hx^2,1 );
    e2DS1=D22+D21+D11;
    
    D33=entS_DisSD3(rho, n, (8*hy)^2,1 );
    D32=entS_DisS32(rho, n, 64*hx*hy^2,2*2 );
    D31=entS_DisS31(rho, n, 32*(hx^2+2)*hy^2,2 );
    D33=real(D33);
    D32=real(D32);
    D31=real(D31);

    e2DS2=D33+D32 +D31 + D22/hy^2*hx^2  +  D21/(hx*hy^4)*(hx^2*hy+2*hy)*hx*hy+ D11/(hy^4*hx^2)*(hx^2*hy+2*hy)^2 ;  
    e2DS1=real(e2DS1);
    e2DS2=real(e2DS2);
    return e2DS1, e2DS2
end


#according to the new error bound in lemma 5, estimate required trotter step r  
function total_step_entropycheck(A, B, n, hx, hy, t, cp, ϵ)
    ϵ1=ϵ/cp;
    t1=t/cp;
    init_left = 1;
    init_right = 1e4;
    r=0;
  
    ComAB=2*(hx*hy*n+hy*2*(n-1));
    ComBBA=4*hx*hy^2*n+8*hy^2*2*(n-1);
    ComAAB=4*hx^2*hy*n+8*hy*(n-1)+8*hy*hx*2*(n-1)+8*hy*(n-2);
    normA=hx*n+n-1;
    normB=hy*n;
  
    for i in 1:cp
  
      P=trotter_error_entropyPF2Para(A, B,n, hx, hy, t1*(i-1));
      e2DS1=P[1];
      e2DS2=P[2];
  
      a1m=16*hx^2*hy^4*n+16*4*hy^4*2*(n-1)+ e2DS1;
      a2m=16*(hx^2+2)^2*hy^2*(n)+16*4*2*hx^2*hy^2*(n-1)+16*4*hy^2*(n-2)+ e2DS2;
      #Error_PF2(r)=r*(sqrt(a1m*(t1/r)^6/144+(5*(t1/r)^7/448*(2*normB+normA)+(t1/r)^8/120*(normB^2+normA^2/4))*ComBBA^2)+ sqrt( a2m*(t1/r)^6/576+(5*(t1/r)^7/1792*(2*normB+normA)+ (t1/r)^8/480*(normB^2+normA^2/4))*ComAAB^2 )+ (t1/r)^4/48*ComBBA+ (t1/r)^4/192*ComAAB);
      Error_PF2(r)=r*(sqrt(a1m*(t1/r)^6/144)+ sqrt( a2m*(t1/r)^6/576)+ (t1/r)^4/32*ComBBA*normA*2+(t1/r)^4/12*ComBBA*normB*2+ (t1/r)^4/32*ComAAB*normB*2+(t1/r)^4/48*ComAAB*normA*2 );
      
      Step= binary_search(Error_PF2, ϵ1, init_left, init_right, 1)
      rworst=total_theo_step_wor2(n, hx, hy, t1, ϵ1);
      Step=min(Step,rworst);
      r=r+Step;
      println(r)
   end
   return r
  end



#according to the new error bound using directly measuring error terms in Eq23 (Method part), estimate required trotter step r  
  function total_step_entropycheck2(A, B, n, hx, hy, t, cp, ϵ)
    ϵ1=ϵ/cp;
    t1=t/cp;
    init_left = 1;
    init_right = 1e4;
    r=0;
  
    ComAB=2*(hx*hy*n+hy*2*(n-1));
    ComBBA=4*hx*hy^2*n+8*hy^2*2*(n-1);
    ComAAB=4*hx^2*hy*n+8*hy*(n-1)+8*hy*hx*2*(n-1)+8*hy*(n-2);
    normA=hx*n+n-1;
    normB=hy*n;
    H=A+B;
    C= A * B - B * A;
    NC1= C*B- B*C; 
    NC2= C*A- A*C;

  
    for i in 1:cp
       psi1=ZeroInp(n);
       a0=expmv(-1im .* t1*(i-1), H, psi1);
       v1=NC1*a0;
       v2=NC2*a0;
      #a1m= opnorm(Array(v1*v1'),2);
      #a2m= opnorm(Array(v2*v2'),2);
      a1m= norm(Array(v1));
      a2m= norm(Array(v2));
      Error_PF2(r)=r*(a1m*(t1/r)^3/12+ a2m*(t1/r)^3/24 + (t1/r)^4/32*ComBBA*normA*2+(t1/r)^4/12*ComBBA*normB*2+ (t1/r)^4/32*ComAAB*normB*2+(t1/r)^4/48*ComAAB*normA*2 );
      
      Step= binary_search(Error_PF2, ϵ1, init_left, init_right, 1)
      rworst=total_theo_step_wor2(n, hx, hy, t1, ϵ1);
      Step=min(Step,rworst);
      r=r+Step;
      println(r)
   end
   return r
  end


function total_theo_step_wor2(n, hx, hy, t, ϵ)

    worst2(r)=r*((4*hx*hy^2*n+8*hy^2*2*(n-1))*(t/r)^3/12 + (4*hx^2*hy*(n)+4*hy*n*2+8*hx*hy*(n-1)*2+8*hy*(n-2))*(t/r)^3/24);
    init_left = 1;
    init_right = 1e4;
    rmin = binary_search( worst2, ϵ, init_left, init_right, 1)
end

function total_theo_step_ave2(n, hx, hy, t, ϵ)
    ave2(r)=r*(sqrt(16*hx^2*hy^4*n+16*4*hy^4*2*(n-1))*(t/r)^3/12+sqrt(16*(hx^2+2)^2*hy^2*(n)+16*4*2*hx^2*hy^2*(n-1)+16*4*hy^2*(n-2))*(t/r)^3/24);
    init_left = 1;
    init_right = 1e4;
    rmin = binary_search( ave2, ϵ, init_left, init_right, 1)
end



function search_empir_num(A, B, t, ϵ)
    #empi_error(r) = empirical_error(A, B, t, r);
#    empi_error(r) = trottervec(A, B, t, r, psi1);
empi_error(r) =  trotter_errorwor2(A, B, t, r);
    init_left = 1;
    init_right = 1e4;
    rmin = binary_search(empi_error, ϵ, init_left, init_right, 1)
end

function search_empir_numRan(A, B, t, ϵ)
    #empi_error(r) = empirical_error(A, B, t, r);
    empi_error(r) = trotter_errorave2(A, B, t, r);
    init_left = 1;
    init_right = 1e4;
    rmin = binary_search(empi_error, ϵ, init_left, init_right, 1)
end


function trotter_errorwor2(A, B, t, r)
    d, ~ = size(H);
    n = log2(d);
    state_ext = zeros(Complex{Float64}, dim);
    state_trot = zeros(Complex{Float64}, dim);
    H=A+B;
    state_ext= exp(-1im .* t*H);
    state_trot=(exp(-1im .* t/r*A/2)*exp(-1im .* t/r*B)*exp(-1im .* t/r*A/2))^r;
    return opnorm(state_ext - state_trot);
end


function trotter_errorave2(A, B, t, r)
    dim, ~ = size(A);
    state_ext = zeros(Complex{Float64}, dim);
    state_trot = zeros(Complex{Float64}, dim);
    gap=zeros(Complex{Float64}, dim);
    H=A+B;
    state_ext= exp(-1im .* t*H);
    state_trot=(exp(-1im .* t/r*A/2)*exp(-1im .* t/r*B)*exp(-1im .* t/r*A/2))^r;
    gap=state_ext - state_trot;
    return  sqrt(tr(gap*gap')/d);
end


 





function commutator!(storage::AbstractSparseArray, A::AbstractSparseArray, B::AbstractSparseArray)
    storage .= A * B - B * A;
end

function binary_search(func, epsilon, init_left, init_right, step)
    # fun is a monotone decreasing function
   left = init_left;
   right = init_right;

   while func(right) > epsilon
       println("Initial guess is too small! 10 times the initial guess!");
       right *= 10;
   end

   if func(left) <= epsilon
       rmin = left;
   else
       while left < right - step
           mid = ceil((left + right)/2);
           if func(mid) > epsilon
               left = mid;
           else
               right = mid;
           end
        #   println("left: ", left, " right: ", right);
       end

       rmin = right;
   end

   return rmin
end

function triangle_errorSe(A::AbstractSparseArray, B::AbstractSparseArray, t, r)
    H = A + B;
    d, ~ = size(H);
    n = log2(d);
    dt = t/r;

    com_AB = spzeros(ComplexF64, d, d);
    com_ABB= spzeros(ComplexF64, d, d);
    com_AAB= spzeros(ComplexF64, d, d);
    commutator!(com_AB, A, B)
    commutator!(com_ABB, com_AB, B)
    commutator!(com_AAB, com_AB, A)

    term1 = sp_trace(com_ABB * com_ABB');
    term1=real(term1);
    term2 = sp_trace(com_AAB * com_AAB');
    term2=real(term2);
    step_err = dt^3 * sqrt(term1/d)/12+ dt^3 * sqrt(term2/d)/24;

    return r * step_err
end



function triangle_errorSeWor(A::AbstractSparseArray, B::AbstractSparseArray, t, r)
    H = A + B;
    d, ~ = size(H);
    n = log2(d);
    dt = t/r;
    com_AB = spzeros(ComplexF64, d, d);
    com_ABB= spzeros(ComplexF64, d, d);
    com_AAB= spzeros(ComplexF64, d, d);
    commutator!(com_AB, A, B)
    commutator!(com_ABB, com_AB, B)
    commutator!(com_AAB, com_AB, A)
  

    term1 = opnorm(Array(com_ABB),2);
    term2 = opnorm(Array(com_AAB),2);
    step_err = dt^3 * term1/12+ dt^3 *term2/24;

    return r * step_err
end


function search_trotter_num_triangle(A, B, t, ϵ)
    #trig_error(r) = triangle_error(A, B, t, r);  #PF1
    trig_error(r) = triangle_errorSe(A, B, t, r); #PF2

    init_left = 1;
    init_right = 1e3;
    rmin = binary_search(trig_error, ϵ, init_left, init_right,1 )
end

function search_trotter_num_triangleWor(A, B, t, ϵ)
    #trig_error(r) = triangle_error(A, B, t, r);  #PF1
    trig_error(r) = triangle_errorSeWor(A, B, t, r); #PF2

    init_left = 1;
    init_right = 1e3;
    rmin = binary_search(trig_error, ϵ, init_left, init_right,1 )
end


function sp_trace(A::AbstractSparseArray)
    m,n = size(A);
    if m != n
        println("Input must be a square matrix. Operation Aborted.");
        return nothing
    else
        tr = 0;
        for i in 1:n
            tr += A[i,i]
        end
        return tr
    end
end

function init_QIMF!(A::AbstractSparseArray, B::AbstractSparseArray, n::Int64, hx::Float64 ,hy::Float64, J::Float64)
    
    # Preparing building blocks
    XX = sparse([0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0]);
    
    X= sparse([0 1;1 0]);
    Y= sparse([0 -1im;1im 0]);
    Z = sparse([1 0;0 -1]);
   

    # Building A and B
    for j in 1:n-1
            left = sparse(Matrix(1*I, 2^(j-1), 2^(j-1)));
            right1 = sparse(Matrix(1*I, 2^(n-j), 2^(n-j)));
            right2 = sparse(Matrix(1*I, 2^(n-1-j), 2^(n-1-j)));

        A .+= hx.*kron(left, X, right1);
        A .+= J.* kron(left, XX, right2);
        B .+= hy.*kron(left, Y, right1);
    end
    left1 = sparse(Matrix(1*I, 2^(n-1), 2^(n-1)));
    A .+= hx.*kron(left1, X);
    B .+= hy.*kron(left1, Y);

end

function von_neumann_entropy(rho::Matrix{Complex{Float64}})
    # 计算密度矩阵的特征值
    eigenvalues = eigvals(rho)

    # 计算冯诺依曼熵
    entropy = 0.0
    for lambda in eigenvalues
        if lambda > 0
            entropy -= lambda * log2(lambda)
        end
    end

    return entropy
end

function TrX4R(rho::Matrix{ComplexF64}, a1::Int64, a2::Int64, a3::Int64, a4::Int64, d::Int64)
   
    G1 = 2^max((a2 - a1 - 1), 0)
    G2 = 2^max(0, (a3 - a2 - 1))
    G3 = 2^max(0, (a4 - a3 - 1))

    M0 = I(2)
    rd = 1
  
    if a2 > a1
        M1 = I(2)
        rd += 1
    else
        M1 = 1
    end

    if a3 > a2
        M2 = I(2)
        rd += 1
    else
        M2 = 1
    end

    if a4 > a3
        M3 = I(2)
        rd += 1
    else
        M3 = 1
    end
    
    M = zeros(Complex{Float64},  2^rd, 2^rd)
   
    for i in 1:2^(a1 - 1)
        for j in 1:G1
            for k in 1:G2
                for p in 1:G3
                    for m in 1:2^(d - a4)
                        A1 = zeros(Complex{Float64}, 1, 2^(a1 - 1))
                        B1 = zeros(Complex{Float64},1, G1)
                        C1 = zeros(Complex{Float64},1, G2)
                        D1 = zeros(Complex{Float64},1, G3)
                        E1 = zeros(Complex{Float64},1, 2^(d - a4))

                        A1[i] = 1
                        B1[j] = 1
                        C1[k] = 1
                        D1[p] = 1
                        E1[m] = 1
                        #println(7)
                        #println(size(V))
                        V = kron(A1, M0, B1, M1, C1, M2, D1, M3, E1)
                        #println(size(V))
                        #println(size(V'))
                        #println(size(V * rho * V' ))
                        #println(size(M))
                        M += V * rho * V' 
                        #println(8)
                    end
                end
            end
        end
    end

    return M
end
   



function TrX6R(rho::Matrix{ComplexF64}, a1::Int64, a2::Int64, a3::Int64, a4::Int64,  a5::Int64,  a6::Int64, d::Int64)
   
    G1 = 2^max((a2 - a1 - 1), 0)
    G2 = 2^max(0, (a3 - a2 - 1))
    G3 = 2^max(0, (a4 - a3 - 1))
    G4 = 2^max(0, (a5 - a4 - 1))
    G5 = 2^max(0, (a6 - a5 - 1))

    M0 = I(2)
    rd = 1
  
    if a2 > a1
        M1 = I(2)
        rd += 1
    else
        M1 = 1
    end

    if a3 > a2
        M2 = I(2)
        rd += 1
    else
        M2 = 1
    end

    if a4 > a3
        M3 = I(2)
        rd += 1
    else
        M3 = 1
    end

    if a5 > a4
        M4 = I(2)
        rd += 1
    else
        M4 = 1
    end

    if a6 > a5
        M5 = I(2)
        rd += 1
    else
        M5 = 1
    end


    
    M = zeros(Complex{Float64},  2^rd, 2^rd)
   
    for i in 1:2^(a1 - 1)
        for j in 1:G1
            for k in 1:G2
                for p in 1:G3
                    for m in 1:G4
                        for ad0 in 1:G5
                            for ad1 in 1:2^(d - a6)
                        A1 = zeros(Complex{Float64}, 1, 2^(a1 - 1))
                        B1 = zeros(Complex{Float64},1, G1)
                        C1 = zeros(Complex{Float64},1, G2)
                        D1 = zeros(Complex{Float64},1, G3)
                        E1 = zeros(Complex{Float64},1, G4)
                        F1=  zeros(Complex{Float64},1, G5)
                        H1=  zeros(Complex{Float64},1, 2^(d - a6))

                        A1[i] = 1
                        B1[j] = 1
                        C1[k] = 1
                        D1[p] = 1
                        E1[m] = 1
                        #println(7)
                        #println(size(V))
                        V = kron(A1, M0, B1, M1, C1, M2, D1, M3, E1)
                        F1[ad0] = 1
                        H1[ad1] = 1
                        V = kron(V, M4,F1, M5, H1)
                        #println(size(V))
                        #println(size(V'))
                        #println(size(V * rho * V' ))
                        #println(size(M))
                        M += V * rho * V' 
                        #println(8)
                            end
                        end
                    end
                end
            end
        end
    end

    return M
end
   



function trace_sparse_complex_matrix(sparse_complex_matrix::SparseMatrixCSC{Complex{Float64}, Int})
    trace_value = Complex{Float64}(0, 0)
    # 遍历对角线元素并累加它们的值
    for i in 1:min(size(sparse_complex_matrix)...)
        trace_value += sparse_complex_matrix[i, i]
    end
    return trace_value
end

function RandInp(n, style)  #1 local random #2 local X #0, haar
    X=[0 1;1 0]
    d=2^n;
  if style==1
  d1=CUE(2)
  Inp=1;
  for j in 1:n
    U1=rand(d1);
    psi=U1*ket(1,2)
   # println(psi)
    Inp=kron(Inp,psi)
  end
  elseif style==2
    Inp=1;
    for j in 1:n
      a=rand(0:1);
      psi=[1; 0]
      psi=X^a*psi;
      #println(psi)
      Inp=kron(Inp,psi)
     end
  else
     se=HaarKet{2}(d);
     Inp=rand(se);
  end
  return Inp
  end
  
  

  function ZeroInp(n)  #1 local random #2 local X #0, haar
   Inp=1;
    for j in 1:n
      psi=[1; 0]
      Inp=kron(Inp,psi)
     end
  
  return Inp
  end
  


  function entS_DisSD3(rho, n, Ex, num)
    entS = 0

    for a1 in 1:(n-2)
        a2 = a1 + 1
        a3=  a1+2 
        for a4 in 1:(n-2)
            a5 = a4 + 1
            a6=  a4 + 2
            b = sort([a1, a2, a3, a4, a5, a6])
            
            rho4 = TrX6R(rho, b[1], b[2], b[3], b[4], b[5], b[6], n)
            d, _ = size(rho4)
            Dis = rho4 - I(d) / d
            ent = tr(sqrt(Dis * Dis'))
            entS = entS + ent
        end
    end
    entS = entS * Ex * num
    return entS
end


function entS_DisS32(rho, n, Ex, num)
    entS = 0

    for a1 in 1:(n-2)
        a2 = a1 + 1
        a3=  a1+2 
        for a4 in 1:(n-1)
            a5 = a4 + 1
            a6=  a5
            b = sort([a1, a2, a3, a4, a5, a6])
            rho4 = TrX6R(rho, b[1], b[2], b[3], b[4], b[5], b[6], n)
            d, _ = size(rho4)
            Dis = rho4 - I(d) / d
            ent = tr(sqrt(Dis * Dis'))
            entS = entS + ent
        end
    end
    entS = entS * Ex * num
    return entS
end




function entS_DisS31(rho, n, Ex, num)
    entS = 0

    for a1 in 1:(n-2)
        a2 = a1 + 1
        a3=  a1+2 
        for a4 in 1:(n)
    
            b = sort([a1, a2, a3, a4])
            rho4 = TrX4R(rho, b[1], b[2], b[3], b[4], n)
            d, _ = size(rho4)
            Dis = rho4 - I(d) / d
            ent = tr(sqrt(Dis * Dis'))
            entS = entS + ent
        end
    end
    entS = entS * Ex * num
    return entS
end


function entS_DisSD2(rho, n, Ex, num)
    entS = 0

    for a1 in 1:(n-1)
        a2 = a1 + 1
        for a3 in 1:(n-1)
            a4 = a3 + 1
            b = sort([a1, a2, a3, a4])
            rho4 = TrX4R(rho, b[1], b[2], b[3], b[4], n)
            d, _ = size(rho4)
            Dis = rho4 - I(d) / d
            ent = tr(sqrt(Dis * Dis'))
            entS = entS + ent
        end
    end
    entS = entS * Ex * num
    return entS
end

function entS_DisS21(rho, n, Ex, num)
    entS = 0

    for a1 in 1:(n-1)
        a2 = a1 + 1
        for a3 in 1:n
            a4 = a3 
            b = sort([a1, a2, a3, a4])
            rho4 = TrX4R(rho, b[1], b[2], b[3], b[4], n)
            d, _ = size(rho4)
            Dis = rho4 - I(d) / d
            ent = tr(sqrt(Dis * Dis'))
            entS = entS + ent
        end
    end
    entS = entS * Ex * num
    return entS
end

 
function entS_DisS11(rho, n, Ex, num)
    entS = 0

    for a1 in 1:n
        a2 = a1 
        for a3 in 1:n
            a4 = a3 
            b = sort([a1, a2, a3, a4])
            rho4 = TrX4R(rho, b[1], b[2], b[3], b[4], n)
            d, _ = size(rho4)
            Dis = rho4 - I(d) / d
            ent = tr(sqrt(Dis * Dis'))
            entS = entS + ent
        end
    end
    entS = entS * Ex * num
    return entS
end


end