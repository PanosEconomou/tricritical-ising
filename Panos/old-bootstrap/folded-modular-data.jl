# Get the modular data for folded minimal models
using LinearAlgebra
using Plots
using PrecompileTools


function FoldedMinimalModularData(p::Int=5,pp::Int=4)

    # Modular data for the minimal model
    cv          = 1 - 6 * (p - pp)^2/ (p * pp)
    hv(r, s)    = ((p*r - pp*s)^2 - (p - pp)^2)/(4*p*pp)
    sv(a, b)    = 2*sqrt(2/(p*pp)) * (isodd(1 + a[2]*b[1] + a[1]*b[2]) ? -1 : 1) * sin(π * p/pp  * a[1]*b[1]) * sin(π * pp/p  * a[2]*b[2])
    tv(a)       = exp(2*π*im*hv(a[1],a[2]) - cv/24)
    labelsv     = vcat(
                       [(r,s) for r in 1:((pp-1)÷2) for s in 1:(p-1)],
                       iseven(pp) ? [(pp÷2,s) for s in 1:((p-1)÷2)] : Tuple{Int,Int}[]
                      )
    index       = Dict(h => i for (i,h) in enumerate(labelsv))
    nv          = length(labelsv)

    # Calculates the labels of the folded primaries
    labels  = vcat(
                  [(labelsv[k÷2],labelsv[k÷2],(-1)^k) for k in 2:(2*nv+1)],
                  [(labelsv[i],labelsv[j],0) for i in 1:nv for j in (i+1):nv],
                  [(labelsv[k÷2],labelsv[k÷2],2*(-1)^k) for k in 2:(2*nv+1)]
                 )

    # Conformal Weights
    function ho(a)
        a1, a2, a3 = a
        δ1 = Int(a == ((1,1),(1,1),-1))
        δ2 = Int(a == ((1,1),(1,1),-2))
        δ3 = Int(abs(a3) == 2)

        return (hv(a1...) + hv(a2...) + a3*(a3 - 1)/2 + δ1)*(1-δ3) + (hv(a1...)/2 + cv/24*(3/2) - (a3 - 2)/8 + δ2)*δ3
    end

    # T-matrix of folded model
    t       = Diagonal([exp(2π*im*(ho(i) - 2*cv/24)) for i in labels])

    # This is the ugliest way I can think of writing the S-matrix 
    svMat   = [sv(a,b) for a in labelsv, b in labelsv]
    tv2     = Diagonal([tv(a)^2 for a in labelsv])
    twist   = svMat*tv2*svMat

    function sMaker(a, b)
        a1, a2, a3 = a
        b1, b2, b3 = b
        aa, ab = abs(a3), abs(b3)
        ia1, ia2 = index[a1], index[a2]
        ib1, ib2 = index[b1], index[b2]
        if     a3 == 0 && b3 == 0;  svMat[ia1,ib1]*svMat[ia2,ib2] + svMat[ia1,ib2]*svMat[ia2,ib1]
        elseif a3 == 0 && ab == 1;  svMat[ia1,ib1]*svMat[ia2,ib1]
        elseif aa == 1 && b3 == 0;  svMat[ia1,ib1]*svMat[ia1,ib2]
        elseif aa == 1 && ab == 1;  svMat[ia1,ib1]^2 / 2
        elseif aa == 1 && ab == 2;  a3 * svMat[ia1,ib1] / 2
        elseif aa == 2 && ab == 1;  b3 * svMat[ia1,ib1] / 2
        elseif aa == 2 && ab == 2
            a3*b3/8 * twist[ia1,ib1] * exp(π*im*(hv(a1...) + hv(b1...) - cv/12))
        else;  zero(ComplexF64)
        end
    end

    s=[sMaker(a,b) for a in labels, b in labels]

    return s, t, labels, svMat
end

using PrecompileTools

@compile_workload begin
    FoldedMinimalModularData(5, 4)
end
