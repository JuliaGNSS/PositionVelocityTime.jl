[![pipeline status](https://git.rwth-aachen.de/nav/PVT.jl/badges/master/pipeline.svg)](https://git.rwth-aachen.de/nav/PVT.jl/commits/master)
[![coverage report](https://git.rwth-aachen.de/nav/PVT.jl/badges/master/coverage.svg)](https://git.rwth-aachen.de/nav/PVT.jl/commits/master)

# Decode GNSS signals.

# Usage
From GNSSDecoder.jl to PVT.jl:

We need decoding results from at least 4 satellites (e.g. decode1, decode2, decode3, decode4)

decRes is an array containing decoding results for navigation data (i.e. decRes = [decode1, decode2, decode3, decode4])

tRXref is the user time when the first preamble was found (not sure if this parameter will be needed)

code_phase is the code phase from tracking results
