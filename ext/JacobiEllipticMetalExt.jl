module JacobiEllipticMetalExt

using JacobiElliptic
using Metal
import JacobiElliptic: CarlsonAlg, FukushimaAlg

# Resolve the process-wide default while constructing a broadcast on the host. The
# resulting concrete callable can then be compiled into a Metal kernel.
for name in JacobiElliptic.func_syms
    f = getfield(JacobiElliptic, name)
    @eval function Base.Broadcast.broadcasted(
        style::Metal.MtlArrayStyle,
        ::typeof($f),
        args...,
    )
        if !isempty(args) && args[1] isa Base.RefValue{<:JacobiElliptic.AbstractAlgorithm}
            return Base.Broadcast.Broadcasted(style, $f, args)
        end

        algorithm =
            JacobiElliptic._effective_default($f, JacobiElliptic.default_algorithm())
        call = JacobiElliptic.AlgorithmCall($f, algorithm)
        return Base.Broadcast.broadcasted(style, call, args...)
    end
end

end
