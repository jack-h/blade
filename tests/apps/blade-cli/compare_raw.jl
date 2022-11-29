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

function compare(i_data, o_data, rtol=0.01)::Result
  correct_count = 0
  count = 0
  if size(i_data) != size(o_data)
    return Result(@sprintf("Shape mismatch: %s != %s", size(i_data), size(o_data)), false)
  end
  dims_correct = Array{Bool}(undef, size(i_data)[2:end])
  for i in CartesianIndices(dims_correct)
    dims_correct[i] = all(isapprox.(real(i_data[:, i]), real(o_data[:, i]), rtol=rtol)) && all(isapprox.(imag(i_data[:, i]), imag(o_data[:, i]), rtol=rtol))
    if !dims_correct[i]
      if count - correct_count < 100
        println(@sprintf("Polarization data mismatch @ %s: %s != %s\n\t(diff: %s)", i, i_data[:, i], o_data[:, i], i_data[:, i] - o_data[:, i]))
      elseif count - correct_count == 100
        println("... omitting more mismatch printouts as 100 have been displayed.")
      end
    else
      correct_count += 1
    end
    count += 1
  end

  Result(@sprintf("%03.06f%% correct (%d/%d)", correct_count/count*100, correct_count, count), all(dims_correct))
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

  rtol = 0.01

  println("\n", compare(i_data, o_data, rtol))

close(i_fio)
close(o_fio)
