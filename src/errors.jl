struct SmallData <: Exception #Data contains not enough Satellites
    msg::String
end


struct InvalidData <: Exception # Size of position array wrong
    msg::String
end


#Decoder errored or warned during decoding or did not decoded properly, position unsafe
struct BadData <: Exception     
    msg::String
end