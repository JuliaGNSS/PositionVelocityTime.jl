struct IncompatibleData <: Exception #Different data array sizes
    msg::String
end


struct SmallData <: Exception #Data contains not enough Satellites
    msg::String
end


struct InvalidData <: Exception # Size of position array wrong
    msg::String
end

struct BrokenData <: Exception #Decoder errored or warned during decoding, position unsafe
    msg::String
end