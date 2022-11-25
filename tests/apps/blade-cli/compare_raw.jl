using Printf
using FFTW
using Blio: GuppiRaw # https://github.com/MydonSolutions/Blio.jl

struct Result
  message::Union{String, Nothing}
  value::Bool
end

function mapToFloat(value::Integer, type::Type)
  return value < 0 ? -1.0*(value / typemin(type)) : value / typemax(type)
end

function mapToFloat(value::Complex{<:Integer}, type::Type)
  return complex(mapToFloat(real(value), real(type)), mapToFloat(imag(value), real(type)))
end

function compare(i_data, o_data, atol=0.01)::Result
  if size(i_data) != size(o_data)
    return Result(@sprintf("Shape mismatch: %s != %s", size(i_data), size(o_data)), false)
  end
  dims_correct = Array{Bool}(undef, size(i_data)[2:end])
  for i in CartesianIndices(dims_correct)
    dims_correct[i] = all(isapprox.(real(i_data[:, i]), real(o_data[:, i]), atol=atol)) && all(isapprox.(imag(i_data[:, i]), imag(o_data[:, i]), atol=atol))
    if !dims_correct[i]
      println(Result(@sprintf("Pol data mismatch @ %s: %s != %s\n\t(atol: %s)", i, i_data[:, i], o_data[:, i], i_data[:, i] - o_data[:, i]), false))
    end
  end

  Result(nothing, all(dims_correct))
end

i_grheader = GuppiRaw.Header()
o_grheader = GuppiRaw.Header()

i_fio = open(ARGS[1], "r")
o_fio = open(ARGS[2], "r")

  read!(i_fio, i_grheader)
  i_data = Array(i_grheader)
  read!(i_fio, i_data)
	if eltype(i_data) <: Complex{<:Integer}
  	i_data = mapToFloat.(i_data, eltype(i_data))
	end

  read!(o_fio, o_grheader)
  o_data = Array(o_grheader)
  read!(o_fio, o_data)

  atol = 0.015

  println(compare(i_data, o_data, atol))

close(i_fio)
close(o_fio)
