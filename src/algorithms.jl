struct SymmetrizedSplitStep <: AbstractNLSEAlgorithm end

abstract type NLSECache <: DECache end
abstract type NLSEConstantCache <: NLSECache end
abstract type NLSEMutableCache <: NLSECache end

mutable struct SplitStepConstantCache <: NLSEConstantCache
    planned_fft!::FFTW.cFFTWPlan
    planned_ifft!::Base.DFT.ScaledPlan
end

function alg_cache(alg::SymmetrizedSplitStep, u, rate_prototype,
                    uEltypeNoUnits, tTypeNoUnits, uprev, uprev2,
                    f, t, dt, reltol, ::Type{Val{false}})
    planned_fft! = plan_fft!(u; flags=FFTW.MEASURE)
    planned_ifft! = plan_ifft!(u; flags=FFTW.MEASURE)
    SplitStepConstantCache(planned_fft!, planned_ifft!)
end
