#include "blade/pipelines/ata/mode_b.hh"

namespace Blade::Pipelines::ATA {

template<typename OT>
ModeB<OT>::ModeB(const Config& config) : config(config), frameJulianDate(1), frameDut1(1) {
    BL_DEBUG("Initializing ATA Pipeline Mode B.");

    if ((config.outputMemPad % sizeof(OT)) != 0) {
        BL_FATAL("The outputMemPad must be a multiple of the output type bytes.")
        BL_CHECK_THROW(Result::ASSERTION_ERROR);
    }

    outputMemPitch = config.outputMemPad + config.outputMemWidth;

    BL_DEBUG("Instantiating input cast from I8 to CF32.");
    this->connect(inputCast, {
        .inputSize = config.numberOfAntennas *
                     config.numberOfFrequencyChannels *
                     config.numberOfTimeSamples *
                     config.numberOfPolarizations,
        .blockSize = config.castBlockSize,
    }, {
        .buf = input,
    });

    BL_DEBUG("Instantiating pre-channelizer with rate {}.", config.preChannelizerRate);
    this->connect(channelizer, {
        .numberOfBeams = 1,
        .numberOfAntennas = config.numberOfAntennas,
        .numberOfFrequencyChannels = config.numberOfFrequencyChannels,
        .numberOfTimeSamples = config.numberOfTimeSamples,
        .numberOfPolarizations = config.numberOfPolarizations,
        .rate = config.preChannelizerRate,
        .blockSize = config.channelizerBlockSize,
    }, {
        .buf = inputCast->getOutput(),
    });

    BL_DEBUG("Instantiating phasor module.");
    this->connect(phasor, {
        .numberOfBeams = config.beamformerBeams,
        .numberOfAntennas = config.numberOfAntennas,
        .numberOfFrequencyChannels = config.numberOfFrequencyChannels * config.preChannelizerRate,
        .numberOfPolarizations = config.numberOfPolarizations,

        .rfFrequencyHz = config.rfFrequencyHz,
        .channelBandwidthHz = config.channelBandwidthHz,
        .totalBandwidthHz = config.totalBandwidthHz,
        .frequencyStartIndex = config.frequencyStartIndex,
        .referenceAntennaIndex = config.referenceAntennaIndex,
        .arrayReferencePosition = config.arrayReferencePosition,
        .boresightCoordinate = config.boresightCoordinate,

        .antennaPositions = config.antennaPositions,
        .antennaCalibrations = config.antennaCalibrations,
        .beamCoordinates = config.beamCoordinates,

        .blockSize = config.phasorsBlockSize,
    }, {
        .frameJulianDate = this->frameJulianDate,
        .frameDut1 = this->frameDut1,
    });

    BL_DEBUG("Instantiating beamformer module.");
    this->connect(beamformer, {
        .numberOfBeams = config.beamformerBeams, 
        .numberOfAntennas = config.numberOfAntennas,
        .numberOfFrequencyChannels = config.numberOfFrequencyChannels * config.preChannelizerRate,
        .numberOfTimeSamples = config.numberOfTimeSamples / config.preChannelizerRate,
        .numberOfPolarizations = config.numberOfPolarizations,
        .enableIncoherentBeam = config.enableIncoherentBeam,
        .enableIncoherentBeamSqrt = false,
        .blockSize = config.beamformerBlockSize,
    }, {
        .buf = channelizer->getOutput(),
        .phasors = phasor->getPhasors(),
    });

    if constexpr (!std::is_same<OT, CF32>::value) {
        BL_DEBUG("Instantiating output cast from CF32 to {}.", typeid(OT).name());
        this->connect(outputCast, {
            .inputSize = beamformer->getOutputSize(),
            .blockSize = config.castBlockSize,
        }, {
            .buf = beamformer->getOutput(),
        });
    }
}

template<typename OT>
Result ModeB<OT>::run(const F64& frameJulianDate,
                      const F64& frameDut1,
                      const Vector<Device::CPU, CI8>& input,
                            Vector<Device::CPU, OT>& output) {
    this->frameJulianDate[0] = frameJulianDate;
    this->frameDut1[0] = frameDut1;

    if (this->getStepCount() == 0) {
        BL_DEBUG("Frame Julian Date: {}", frameJulianDate);
        BL_DEBUG("Frame DUT1: {}", frameDut1);
    }

    BL_CHECK(this->copy(inputCast->getInput(), input));
    BL_CHECK(this->compute());
    BL_CHECK(this->copy2D(
        output,
        outputMemPitch,         // dpitch
        0,
        this->getOutput(),      // src
        config.outputMemWidth,  // spitch
        0,
        config.outputMemWidth,  // width
        (beamformer->getOutputSize() * sizeof(OT)) / config.outputMemWidth
    ));

    return Result::SUCCESS;
}

template<typename OT>
Result ModeB<OT>::run(const F64& frameJulianDate,
                      const F64& frameDut1,
                      const Vector<Device::CPU, CI8>& input, 
                      const U64& outputBlockIndex,
                      const U64& outputNumberOfBlocks,
                            Vector<Device::CUDA, OT>& output) {
    this->frameJulianDate[0] = frameJulianDate;
    this->frameDut1[0] = frameDut1;

    if (this->getStepCount() == 0) {
        BL_DEBUG("Frame Julian Date: {}", frameJulianDate);
        BL_DEBUG("Frame DUT1: {}", frameDut1);
    }

    BL_CHECK(this->copy(inputCast->getInput(), input));
    BL_CHECK(this->compute());

    const auto& width = (beamformer->getOutputSize() / config.beamformerBeams / 
            (config.numberOfFrequencyChannels * config.preChannelizerRate)) * sizeof(OT);
    const auto& height = config.beamformerBeams * 
            (config.numberOfFrequencyChannels * config.preChannelizerRate);

    BL_CHECK(this->copy2D(
        output,
        width * outputNumberOfBlocks,
        0,
        this->getOutput(),
        width,
        0,
        width,
        height
    )); 

    return Result::SUCCESS;
}

template class BLADE_API ModeB<CF16>;
template class BLADE_API ModeB<CF32>;

}  // namespace Blade::Pipelines::ATA
