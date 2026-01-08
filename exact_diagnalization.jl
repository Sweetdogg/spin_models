# 对角化张量积形式的 TFIM 哈密顿量
using Printf, LinearAlgebra, LinearMaps, OMEinsum, KrylovKit
function doApplyHam(psiIn::AbstractVector, hloc, N)
    d = size(hloc,1)
    psiOut = zeros(eltype(psiIn), d^N)
    # hloc[a,b,i,j] 作用在 (s_k=i, s_{k+1}=j)
    @inbounds for k in 1:N-1
        ψ = reshape(psiIn, d^(k-1), d, d, d^(N-1-k))
        psiOut .+= vec(ein"abij,lijr->labr"(hloc, ψ))
    end
    ψ = reshape(psiIn, d, d^(N-2), d)            # (s1, mid, sN)
    psiOut .+= vec(ein"abij,jmi->bma"(hloc, ψ))   # PBC
    return psiOut
end

N = 20 # 尺寸
numval = 2 # 计算最低能量本征值个数
d = 2 # 局域维度
dim = d^N
sX = [0 1; 1 0]; sZ = [1 0; 0 -1]; sI = [1 0; 0 1]
g = 1.0
hloc = reshape(-kron(sZ,sZ) - 0.5*g*kron(sI,sX) - 0.5*g*kron(sX,sI),2,2,2,2)
EnExact = -2/sin(pi/(2*N))

Hmap = LinearMap(ψ -> doApplyHam(ψ, hloc, N), dim;
    ismutating=false, issymmetric=true, ishermitian=true, isposdef=false)

diagtime = @elapsed begin
    Energy, vecs, info = eigsolve(Hmap, dim, numval, :SR;
                                krylovdim=20,maxiter=300, tol=1e-15)
    psi = hcat(vecs...)
end

Err = Energy[1] - EnExact
Δ = Energy[2] - Energy[1]
@printf "N: %d Time: %.3f Energy: %e Err: %e Gap: %e\n" N diagtime Energy[1] Err Δ

m_2(psi::AbstractMatrix, N::Int; d::Int=2) = begin
    m_z = [(2*count_ones(s) - N)/N for s in 0:d^N-1]
    prob = abs2.(psi[:,1])
    dot(prob, m_z .^ 2)
end
@printf "g=%f  <m²>=%e\n" g m_2(psi,N;d=d)

ent_entropy(ψ::AbstractVector, d::Int, l::Int) = begin
    s = svdvals(reshape(ψ, d^l, :))
    p = s .^ 2
    -sum(p .* log.(p))
end
l = N ÷ 2
S_half = ent_entropy(psi[:,1], d, l)
@printf "Entanglement Entropy S(%d) at l=%d: %e\n" N l S_half