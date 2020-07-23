module PVT

    using DocStringExtensions, Parameters, FixedPointNumbers, LinearAlgebra, GNSSDecoder, GNSSSignals
    using Unitful: s, Hz
    

    #TODO: PVT Interface here
    # inputs: Array of Decoder, Code_phases and eventually Carrier_phases
    #
    #
    #
    #
    #output: Struct of User Position and DOP
    
    include("pseudo_range.jl")
    include("user_position.jl")
    include("sat_position.jl")
    include("sv_time.jl")
    include("errors.jl")
end #end module
