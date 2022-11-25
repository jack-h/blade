import os, sys, re, argparse
import numpy
import pyproj
from guppi.guppi import Guppi # https://github.com/MydonSolutions/guppi/tree/write

import astropy.constants as const
from astropy.coordinates import ITRS, SkyCoord
from astropy.time import Time
import astropy.units as u

def upchannelize_frequencies(frequencies, rate):
    fine_frequencies = numpy.zeros((len(frequencies), rate), dtype=numpy.float64)
    chan_bw = 0
    for coarse_chan_i in range(len(frequencies)-1):
        chan_bw = frequencies[coarse_chan_i+1] - frequencies[coarse_chan_i]
        fine_frequencies[coarse_chan_i, :] = numpy.linspace(
            frequencies[coarse_chan_i],
            frequencies[coarse_chan_i+1],
            rate,
            endpoint=False
        )
    fine_frequencies[-1, :] = numpy.linspace(
        frequencies[-1],
        frequencies[-1]+chan_bw,
        rate,
        endpoint=False
    )
    return fine_frequencies.flatten()


def upchannelize(
    datablock: numpy.ndarray, # [Antenna, Frequency, Time, Polarization]
    rate: int
):
    """
    Params
    ------
        datablock: numpy.ndarray # [Antenna, Frequency, Time, Polarization]
            the data to be upchannelized
        frequencies: list
            the frequency of each channel, will be upchannelized and returned
        rate: int
            the FFT length
    
    Return
    ------
    upchannelized_datablock, upchannelized_frequencies
    """
    # Be Faster!
    A, F, T, P = datablock.shape
    assert T % rate == 0, f"Rate {rate} is not a factor of time {T}."
    fine_datablock = datablock.copy()
    fine_datablock = fine_datablock.reshape((A, F, T//rate, rate, P))
    fine_datablock = numpy.fft.fftshift(numpy.fft.fft(
            fine_datablock,
            axis=3
        ),
        axes=3
    )
    fine_datablock = numpy.transpose(fine_datablock, (0, 1, 3, 2, 4))
    return fine_datablock.reshape((A, F*rate, T//rate, P))

    # not slower!
    A, F, T, P = datablock.shape
    assert T % rate == 0, f"Rate {rate} is not a factor of time {T}."
    output = numpy.zeros((A, F*rate, T//rate, P), dtype=numpy.complex64)
    for iant in range(A):
        ant = datablock[iant]

        for ichan in range(F):
            time_pol = ant[ichan]

            for ipol in range(P):
                time_arr = time_pol[:, ipol]

                for ispec in range(output.shape[2]):
                    output[
                        iant,
                        ichan*rate : (ichan+1)*rate,
                        ispec,
                        ipol
                    ] = numpy.fft.fftshift(numpy.fft.fft(
                        time_arr[ispec*rate:(ispec+1)*rate]
                    ))

    return output


def _compute_uvw(ts, source, ant_coordinates, ref_coordinates):
    """Computes UVW antenna coordinates with respect to reference

    Copyright 2021 Daniel Estevez <daniel@destevez.net>

    Args:
        ts: array of Times to compute the coordinates
        source: source SkyCoord
        ant_coordinates: antenna ECEF coordinates.
            This is indexed as (antenna_number, xyz)
        ref_coordinates: phasing reference centre coordinates.
            This is indexed as (xyz)

    Returns:
        The UVW coordinates in metres of each of the baselines formed
        between each of the antennas and the phasing reference. This
        is indexed as (time, antenna_number, uvw)
    """
    baselines_itrs = ant_coordinates - ref_coordinates

    # Calculation of vector orthogonal to line-of-sight
    # and pointing due north.
    north_radec = [source.ra.deg, source.dec.deg + 90]
    if north_radec[1] > 90:
        north_radec[1] = 180 - north_radec[1]
        north_radec[0] = 180 + north_radec[0]
    north = SkyCoord(ra=north_radec[0]*u.deg, dec=north_radec[1]*u.deg)

    source_itrs = source.transform_to(ITRS(obstime=Time(ts))).cartesian
    north_itrs = north.transform_to(ITRS(obstime=Time(ts))).cartesian
    east_itrs = north_itrs.cross(source_itrs)

    ww = baselines_itrs @ source_itrs.xyz.value
    vv = baselines_itrs @ north_itrs.xyz.value
    uu = baselines_itrs @ east_itrs.xyz.value
    uvw = numpy.stack((uu.T, vv.T, ww.T), axis=-1)

    return uvw


def _create_delay_phasors(delay, frequencies):
    return numpy.exp(-1j*2*numpy.pi*delay*frequencies)


def _get_fringe_rate(delay, fringeFrequency):
    return numpy.exp(-1j*2*numpy.pi*delay * fringeFrequency)


def phasors(
    antennaPositions: numpy.ndarray, # [Antenna, XYZ]
    boresightCoordinate: SkyCoord, # degrees
    beamCoordinates: 'list[SkyCoord]', #  degrees
    times: numpy.ndarray, # [unix]
    frequencies: numpy.ndarray, # [channel-frequencies] Hz
    calibrationCoefficients: numpy.ndarray, # [Frequency-channel, Polarization, Antenna]
    referenceAntennaIndex: int = 0
):
    """
    Return
    ------
        phasors (B, A, F, T, P), delays (T, A, B)

    """

    assert frequencies.shape[0] % calibrationCoefficients.shape[0] == 0, f"Calibration Coefficients' Frequency axis is not a factor of frequencies: {calibrationCoefficients.shape[0]} vs {frequencies.shape[0]}."

    phasorDims = (
        beamCoordinates.shape[0],
        antennaPositions.shape[0],
        frequencies.shape[0],
        times.shape[0],
        calibrationCoefficients.shape[1]
    )
    calibrationCoeffFreqRatio = frequencies.shape[0] // calibrationCoefficients.shape[0]

    phasors = numpy.zeros(phasorDims, dtype=numpy.complex128)
    
    delays = numpy.zeros(
        (
            times.shape[0],
            antennaPositions.shape[0],
            beamCoordinates.shape[0],
        ),
        dtype=numpy.float64
    )

    for t, tval in enumerate(times):
        ts = Time(tval, format='unix')
        boresightUvw = _compute_uvw(
            ts,
            boresightCoordinate, 
            antennaPositions,
            antennaPositions[referenceAntennaIndex],
        )
        for b in range(phasorDims[0]):
            # These UVWs are centred at the reference antenna, 
            # i.e. UVW_irefant = [0, 0, 0]
            beamUvw = _compute_uvw( # [Antenna, UVW]
                ts,
                beamCoordinates[b], 
                antennaPositions,
                antennaPositions[referenceAntennaIndex],
            )

            antennaBeamDelays = (beamUvw[:,2] - boresightUvw[:,2]) / const.c.value
            delays[t, :, b] = antennaBeamDelays
            for a, delay in enumerate(antennaBeamDelays):
                phasors[b, a, :, t, 0] = _create_delay_phasors(
                    delay,
                    frequencies
                )
                phasors[b, a, :, t, 0] *= _get_fringe_rate(
                    delay,
                    frequencies[0]
                )
                for p in range(phasorDims[-1]):
                    # Phasor is the same for each polarization, until the 
                    # calibration coefficent makes it different.
                    # Reuse phasor of pol0, changing it last
                    revP = phasorDims[-1]-1-p
                    for c in range(calibrationCoeffFreqRatio):
                        fine_slice = range(c, frequencies.shape[0], calibrationCoeffFreqRatio)
                        phasors[b, a, fine_slice, t, revP] = phasors[b, a, fine_slice, t, 0] * calibrationCoefficients[:, revP, a]
    return phasors, delays


def beamform(
    datablock: numpy.ndarray, # [Antenna, Frequency, Time, Polarization]
    phasors: numpy.ndarray, # [Beam, Aspect, Frequency, Time=1, Polarization]
    lastBeamIncoherent = False
):
    # Beamform with Numpy.
    output = numpy.multiply(datablock, phasors) 
    
    if lastBeamIncoherent:
        output[-1] = output[-1].real**2 + output[-1].imag**2

    output = output.sum(axis=1) # sum across antenna, collapsing that dimension

    if lastBeamIncoherent:
        output[-1] = numpy.sqrt(output[-1])
    return output


def guppiheader_get_blockdims(guppiheader):
    nof_obschan = guppiheader.get('OBSNCHAN', 1)
    nof_ant = guppiheader.get('NANTS', 1)
    nof_chan = nof_obschan // nof_ant
    nof_pol = guppiheader.get('NPOL', 1)
    nof_time = guppiheader.get('BLOCSIZE', 0) * 8 // (nof_ant * nof_chan * nof_pol * 2 * guppiheader.get('NBITS', 8))
    return (nof_ant, nof_chan, nof_time, nof_pol)

def guppiheader_get_unix_midblock(guppiheader):
    dims = guppiheader_get_blockdims(guppiheader)
    ntime = dims[2]
    synctime = guppiheader.get('SYNCTIME', 0)
    pktidx = guppiheader.get('PKTIDX', 0)
    tbin = guppiheader.get('TBIN', 1.0/guppiheader.get("CHAN_BW", 0.0))
    piperblk = guppiheader.get('PIPERBLK', ntime)
    return synctime + pktidx * ((tbin * ntime)/piperblk)

def index_of_time(times, t):
    for i, ti in enumerate(times):
        if ti == t:
            return i
        if ti > t:
            assert i != 0, f"Time {t} appears to be before the start of times: {times[0]}"
            return i

    assert False, f"Time {t} appears to be past the end of times: {times[-1]}"


if __name__ == "__main__":
    import h5py

    parser = argparse.ArgumentParser(
        description='Generate a (RAW, BFR5):(Filterbank) input:output pair of beamforming',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('bfr5', type=str)
    parser.add_argument('guppi', type=str)

    parser.add_argument('-u', '--upchannelization-rate', type=int, default=1,)
    parser.add_argument('-o', '--output-directory', type=str, default=None,)

    args = parser.parse_args()
    bfr5 = h5py.File(args.bfr5, 'r')
    gfile = Guppi(args.guppi)

    guppi_stem = re.match(r"(.*)\.\d{4}.raw", os.path.basename(args.guppi)).group(1)
    if args.output_directory is None:
        args.output_directory = os.path.dirname(args.guppi)
    
    datablockDims = (
        bfr5['diminfo']['nants'][()],
        bfr5['diminfo']['nchan'][()],
        bfr5['diminfo']['ntimes'][()],
        bfr5['diminfo']['npol'][()],
    )
    
    beamCoordinates = numpy.array([
        SkyCoord(
            bfr5['beaminfo']['ras'][beamIdx],
            bfr5['beaminfo']['decs'][beamIdx],
            unit='rad'
        )
        for beamIdx in range(bfr5['diminfo']['nbeams'][()])
    ])
    boresightCoordinate = SkyCoord(
        bfr5['obsinfo']['phase_center_ra'][()],
        bfr5['obsinfo']['phase_center_dec'][()],
        unit='rad'
    )

    antennaPositions = bfr5['telinfo']['antenna_positions'][:]
    antennaPositionFrame = bfr5['telinfo']['antenna_position_frame'][()]

    # convert from antennaPositionFrame to 'xyz'
    if antennaPositionFrame == 'ecef':
        transformer = pyproj.Proj.from_proj(
            pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84'),
            pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84'),
        )
        lla = (
            bfr5['telinfo']['longitude'][()],
            bfr5['telinfo']['latitude'][()],
            bfr5['telinfo']['altitude'][()],
        )
        telescopeCenterXyz = transformer.transform(*lla,radians=False)
        print(f"Subtracting Telescope center from ECEF antenna-positions: {telescopeCenterXyz}")
        for i in range(antennaPositions.shape[0]):
            antennaPositions[i, :] -= telescopeCenterXyz

    frequencies = bfr5['obsinfo']['freq_array'][:] * 1e9
    times = bfr5['delayinfo']['time_array'][:]

    upchan_frequencies = upchannelize_frequencies(
        frequencies,
        args.upchannelization_rate
    )

    phasorCoeffs, delays = phasors(
        antennaPositions, # [Antenna, XYZ]
        boresightCoordinate, # degrees
        beamCoordinates, #  degrees
        times, # [unix]
        upchan_frequencies, # [channel-frequencies] Hz
        bfr5['calinfo']['cal_all'][:], # [Antenna, Frequency-channel, Polarization]
        referenceAntennaIndex = 0
    )

    recipe_delays_agreeable = bfr5['delayinfo']['delays'][:] == delays # should probably use isclose
    if not recipe_delays_agreeable.all():
        print(f"The delays in the provided recipe file do not match the calculated delays:\n{recipe_delays_agreeable}")
        print(f"Using calculated delays:\n{delays}")
    else:
        print("The recipe file's delays match the calculated delays.")

    hdr, data = gfile.read_next_block()
    block_index = 1
    block_times_index = 0
    file_open_mode = 'wb'
    while True:
        if hdr is None:
            break

        assert datablockDims == guppiheader_get_blockdims(hdr), f"#{block_index}: {datablockDims} != {guppiheader_get_blockdims(hdr)}"

        upchan_data = upchannelize(
            data,
            args.upchannelization_rate
        )

        try:
            block_times_index = index_of_time(
                times,
                guppiheader_get_unix_midblock(hdr)
            )
        except:
            pass
        beams = beamform(
            upchan_data,
            phasorCoeffs[:, :, :, block_times_index:block_times_index+1, :]
        )
        
        hdr["DATATYPE"] = "FLOAT"
        hdr["TBIN"] *= args.upchannelization_rate
        hdr["CHAN_BW"] /= args.upchannelization_rate
        for beam in range(beams.shape[0]):
            Guppi.write_to_file(
                os.path.join(args.output_directory, f"{guppi_stem}-beam{beam:03d}.0000.raw"),
                hdr,
                beams[beam:beam+1, ...].astype(numpy.complex64),
                file_open_mode=file_open_mode
            )
        file_open_mode = 'ab'
        hdr, data = gfile.read_next_block()
        block_index += 1
