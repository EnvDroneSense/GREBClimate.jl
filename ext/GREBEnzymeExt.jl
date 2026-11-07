# Enzyme derivative rules for the fast scalar log/exp used in the plain loops
# of LWradiation! and hydro!. Loaded only when Enzyme is loaded.
module GREBEnzymeExt

using GREBClimate: fastlog, fastexp
using Enzyme.EnzymeRules
using Enzyme: Const, Active, Duplicated, DuplicatedNoNeed

const Dup = Union{Duplicated,DuplicatedNoNeed}

# The primal goes into an augmented return only when Enzyme asks for it
primal_if_needed(config, value) = needs_primal(config) ? value : nothing

# What a forward rule returns depends on whether Enzyme wants the primal, the shadow or both
function forward_return(config, primal, shadow)
    needs_primal(config) && needs_shadow(config) && return Duplicated(primal, shadow)
    needs_shadow(config) && return shadow
    return primal
end

# d log(x) = dx / x
EnzymeRules.forward(config::FwdConfig, ::Const{typeof(fastlog)}, ::Type{<:Dup}, x::Duplicated) = forward_return(config, fastlog(x.val), x.dval / x.val)
EnzymeRules.forward(::FwdConfig, ::Const{typeof(fastlog)}, ::Type{<:Const}, x::Const) = fastlog(x.val)
EnzymeRules.augmented_primal(config::RevConfigWidth{1}, ::Const{typeof(fastlog)}, ::Type{<:Active}, x::Active) = AugmentedReturn(primal_if_needed(config, fastlog(x.val)), nothing, nothing)
EnzymeRules.reverse(::RevConfigWidth{1}, ::Const{typeof(fastlog)}, dret::Active, tape, x::Active) = (dret.val / x.val,)

# d exp(x) = exp(x) dx; the reverse pass reuses exp(x) stored as the tape
EnzymeRules.forward(config::FwdConfig, ::Const{typeof(fastexp)}, ::Type{<:Dup}, x::Duplicated) = (e=fastexp(x.val); forward_return(config, e, x.dval * e))
EnzymeRules.forward(::FwdConfig, ::Const{typeof(fastexp)}, ::Type{<:Const}, x::Const) = fastexp(x.val)
EnzymeRules.augmented_primal(config::RevConfigWidth{1}, ::Const{typeof(fastexp)}, ::Type{<:Active}, x::Active) = (e=fastexp(x.val); AugmentedReturn(primal_if_needed(config, e), nothing, e))
EnzymeRules.reverse(::RevConfigWidth{1}, ::Const{typeof(fastexp)}, dret::Active, tape, x::Active) = (dret.val * tape,)

# Constant arguments (for example a temperature that is not differentiated)
for f in (fastlog, fastexp)
    @eval EnzymeRules.augmented_primal(config::RevConfigWidth{1}, ::Const{typeof($f)}, ::Type{<:Const}, x::Const) = AugmentedReturn(primal_if_needed(config, $f(x.val)), nothing, nothing)
    @eval EnzymeRules.reverse(::RevConfigWidth{1}, ::Const{typeof($f)}, ::Type{<:Const}, tape, x::Const) = (nothing,)
end

end
