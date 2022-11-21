module MINE

using MINE_jll

const EST_MIC_APPROX = Cint(0)
const EST_MIC_E = Cint(1)

EST = Dict(
    "mic_approx" => EST_MIC_APPROX,
    "mic_e" => EST_MIC_E
)

struct MINEProblem
    n::Cint
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
end

MINEProblem(x::Vector, y::Vector) = MINEProblem(length(x), pointer(x), pointer(y))

struct MINEParameter
    alpha::Cdouble
    c::Cdouble
    est::Cint
end

struct MINEScore
    n::Cint  # number of rows of M
    m::Ptr{Cint}  # number of cols of M[i] for each i
    M::Ptr{Ptr{Cdouble}}  # the (equi)characteristic matrix
end

struct MINEMatrix
    data::Ptr{Cdouble}  # matrix in row-major order
    n::Cint  # number of rows
    m::Cint  # number of cols
end

struct MINEpstats
    mic::Ptr{Cdouble} # condensed matrix
    tic::Ptr{Cdouble} # condensed matrix
    n::Cint  # number of elements
end

struct MINEcstats
    mic::Ptr{Float64}  # matrix in row-major order
    tic::Ptr{Float64}  # matrix in row-major order
    n::Cint  # number of rows
    m::Cint  # number of cols
end


function mine_check_parameter(mine_parameter::MINEParameter)
    return ccall(
        (:mine_check_parameter, libmine),
        Cchar,
        (Ptr{MINEParameter},),
        Ref(mine_parameter),
    )
end

function mine_compute_score(mine_problem::MINEProblem, mine_parameter::MINEParameter)::Ptr{MINEScore}
    return ccall(
        (:mine_compute_score, libmine),
        Ptr{MINEScore},
        (Ptr{MINEProblem}, Ptr{MINEParameter}),
        Ref(mine_problem),
        Ref(mine_parameter),
    )
end

function mine_compute_score(x::Vector, y::Vector; alpha=0.6, c=15, est="mic_approx")
    if length(x) != length(y)
        error("x, y: shape mismatch")
    end
    mine_problem = MINEProblem(x, y)
    mine_parameter = MINEParameter(alpha, c, EST[est])
    return mine_compute_score(mine_problem, mine_parameter)
end

function get_score(mine_score::Ptr{MINEScore})
    M = Vector{Float64}[]
    for i in range(1, unsafe_load(mine_score).n)
        M_i = unsafe_wrap(
            Vector{Float64},
            unsafe_load(unsafe_load(mine_score).M, i),
            unsafe_load(unsafe_load(mine_score).m, i)
        )
        push!(M, M_i)
    end
    return M
end

function mine_mic(mine_score::Ptr{MINEScore})
    return ccall(
        (:mine_mic, libmine),
        Cdouble,
        (Ptr{MINEScore},),
        mine_score,
    )
end

function mine_mas(mine_score::Ptr{MINEScore})
    return ccall(
        (:mine_mas, libmine),
        Cdouble,
        (Ptr{MINEScore},),
        mine_score,
    )
end

function mine_mev(mine_score::Ptr{MINEScore})
    return ccall(
        (:mine_mev, libmine),
        Cdouble,
        (Ptr{MINEScore},),
        mine_score,
    )
end

function mine_mcn(mine_score::Ptr{MINEScore}; eps=0)
    return ccall(
        (:mine_mcn, libmine),
        Cdouble,
        (Ptr{MINEScore}, Cdouble),
        mine_score,
        eps,
    )
end
mine_mcn(mine_score, eps=0) = mine_mcn(mine_score; eps=eps)

function mine_mcn_general(mine_score::Ptr{MINEScore})
    return ccall(
        (:mine_mcn_general, libmine),
        Cdouble,
        (Ptr{MINEScore},),
        mine_score,
    )
end

function mine_tic(mine_score::Ptr{MINEScore}; norm=false)
    return ccall(
        (:mine_tic, libmine),
        Cdouble,
        (Ptr{MINEScore}, Cint),
        mine_score,
        norm
    )
end
mine_tic(mine_score, norm=false) = mine_tic(mine_score; norm=norm)

function mine_gmic(mine_score::Ptr{MINEScore}; p=-1)
    return ccall(
        (:mine_gmic, libmine),
        Cdouble,
        (Ptr{MINEScore}, Cdouble),
        mine_score,
        p
    )
end
mine_gmic(mine_score, p=-1) = mine_gmic(mine_score; p=p)

# TODO
function pstats(X; alpha=0.6, c=15, est="mic_approx")

    param = MINEParameter(alpha, c, EST[est])
    ret = mine_check_parameter(param)
    if ret != 0
       error(ret)
    end
    Xa = convert(Array, X)
    Xm = MINEMatrix(pointer(Xa), size(Xa)[1], size(Xa)[2])

    pstats_ptr = ccall(
        (:mine_compute_pstats, libmine),
        Ptr{MINEpstats},
        (Ptr{MINEMatrix}, Ptr{MINEParameter}),
        Ref(Xm),
        Ref(param),
    )
    pstats = unsafe_load(pstats_ptr)
    mic = unsafe_wrap(
        Vector{Float64},
        pstats.mic,
        pstats.n,
    )
    tic = unsafe_wrap(
        Vector{Float64},
        pstats.tic,
        pstats.n,
    )

    return mic, tic
end

# TODO
function cstats(X, Y; alpha=0.6, c=15, est="mic_approx")

    if size(X)[2] != size(X)[2]
        error("X, Y: shape mismatch")
    end

    param = MINEParameter(alpha, c, EST[est])
    ret = mine_check_parameter(param)
    if ret != 0
       error(ret)
    end

    Xa = convert(Array, X)
    Ya = convert(Array, Y)
    Xm = MINEMatrix(pointer(Xa), size(Xa)[1], size(Xa)[2])
    Ym = MINEMatrix(pointer(Ya), size(Ya)[1], size(Ya)[2])
    cstats_ptr = ccall(
        (:mine_compute_cstats, libmine),
        Ptr{MINEcstats},
        (Ptr{MINEMatrix}, Ptr{MINEMatrix}, Ptr{MINEParameter}),
        Ref(Xm),
        Ref(Ym),
        Ref(param),
    )
    cstats = unsafe_load(cstats_ptr)

    mic = Matrix{Float64}(undef, cstats.n, cstats.m)
    tic = Matrix{Float64}(undef, cstats.n, cstats.m)
    for i in range(1, cstats.n)
        for j in range(1, cstats.m)
            mic[i, j] = unsafe_load(cstats.mic, ((i - 1) * cstats.n) + j)
            tic[i, j] = unsafe_load(cstats.tic, ((i - 1) * cstats.n) + j)
        end
    end

    return mic, tic
end

end  # module MINE
