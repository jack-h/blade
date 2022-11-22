include("helper.jl")
using Dates

using ATA_BFR5_Genie # https://github.com/MydonSolutions/ata_bfr5_genie

# Zeros BFR5, with `1+0j` for cal_all
# Complex-Exponential Signal repeated throughout RAW data

now_unix = datetime2unix(now())

generateTestInputs(
	"synthetic_test_rand",
	function (recipe)
		n_ant = recipe.diminfo.nants
		n_beam = recipe.diminfo.nbeams
		n_chan_perant = recipe.diminfo.nchan
		
		recipe.telinfo.antenna_positions = rand(Float64, (3, recipe.diminfo.nants)) .* 10000.0
		
		recipe.beaminfo.ras = rand(Float64, (recipe.diminfo.nbeams)) .* 0.5
		recipe.beaminfo.decs = rand(Float64, (recipe.diminfo.nbeams)) .* 0.5

		recipe.obsinfo.freq_array = collect(0:recipe.diminfo.nchan-1) .* 0.5e6 .+ 1e6
		recipe.obsinfo.freq_array /= 1e9 # stored in GHz
		recipe.obsinfo.phase_center_ra = recipe.beaminfo.ras[1]
		recipe.obsinfo.phase_center_dec = recipe.beaminfo.decs[1]

		recipe.delayinfo.time_array = zeros(Float64, (1)) .+ now_unix
		recipe.delayinfo.jds = [
			(unix / 86400) + 2440587.5
			for unix in recipe.delayinfo.time_array
		]

		recipe.calinfo.cal_all = rand(ComplexF32, (recipe.diminfo.nants, recipe.diminfo.npol, recipe.diminfo.nchan))

		recipe.delayinfo.delays = cat(
			(
				ATA_BFR5_Genie.calculateBeamDelays(
					recipe.telinfo.antenna_positions,
					1,
					recipe.beaminfo.ras[1], recipe.beaminfo.decs[1],
					hcat(recipe.beaminfo.ras, recipe.beaminfo.decs)',
					recipe.telinfo.longitude, recipe.telinfo.latitude, recipe.telinfo.altitude,
					unix, recipe.delayinfo.dut1
				)
				for unix in recipe.delayinfo.time_array
			)...
			;
			dims=3
		)

	end,
	function (header, data)
		n_ant = header["NANTS"]
		n_chan_perant = div(header["OBSNCHAN"], n_ant)
		n_time = header["PIPERBLK"]
		n_pol = header["NPOL"]

		header["SYNCTIME"] = now_unix

		data .= rand(eltype(data), (n_pol, n_time, n_chan_perant, n_ant))
	end,
	directory = "/mnt/buf1/mydonsol_blade/basics",
	iterateRawcallback = true
)
