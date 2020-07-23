struct IncompatibleData <: Exception
    msg::String
end


struct SmallData <: Exception
    msg::String
end


struct InvalidData <: Exception
    msg::String
end