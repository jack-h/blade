# Fundamental Beamforming Tests

For the purposes of brevity, the RAW data is spoken about programmatically, with dimensions [fastest=pol, time, chan, slowets=ant], as are the antenna coefficients, with dimensions [fastest=ant, pol, slowest=chan].

The synthesized tests have the following RAW block dimension, unless otherwise stated: [A: 4, F: 128, T: 16384, P: 2, 8bit-complex].

One coherent and one incoherent beam are formed, unless otherwise stated;

The BLADE-CLI is run with `-C 64` to ensure that the coarse frequency-channel stepping works, and with `-T 4096`, unless otherwise stated.

Each test is collapsed below.

<details><summary>0. cal_all Ones, delays Zeros, RAW signal in [:, :, NCHAN/2, :]</summary>


<details><summary>GUPPI RAW Input</summary>

![synthetic_test_0](./plots/synthetic_test_0.png)

</details>

<details><summary>Beamformed Output (No upchannelization)</summary>

![synthetic_test_0_c1_beam0](./plots/synthetic_test_0_c1_beam0.png)
</details>

<details><summary>Beamformed Output (upchannelization rate of 4)</summary>

![synthetic_test_0_c4_beam0](./plots/synthetic_test_0_c4_beam0.png)
![synthetic_test_0_c4_beam0_zoom](./plots/synthetic_test_0_c4_beam0_zoom.png)
</details>


</details>


<details><summary>1. cal_all Ones, delays Zeros, RAW signal in [:, :, NCHAN/2, NANT/2]</summary>


<details><summary>GUPPI RAW Input</summary>

![synthetic_test_1](./plots/synthetic_test_1.png)

</details>

<details><summary>Beamformed Output (No upchannelization)</summary>

![synthetic_test_1_c1_beam0](./plots/synthetic_test_1_c1_beam0.png)
</details>

<details><summary>Beamformed Output (upchannelization rate of 4)</summary>

![synthetic_test_1_c4_beam0](./plots/synthetic_test_1_c4_beam0.png)
</details>


</details>


<details><summary>2. delays Zeros, cal_all band pass in [:, :, NCHAN/2], RAW signal in [:, :, NCHAN/2, NANT/2] and [:, :, 3*NCHAN/2, :]</summary>


<details><summary>GUPPI RAW Input</summary>

![synthetic_test_2](./plots/synthetic_test_2.png)

</details>

<details><summary>Beamformed Output (No upchannelization)</summary>

![synthetic_test_2_c1_beam0](./plots/synthetic_test_2_c1_beam0.png)
</details>

<details><summary>Beamformed Output (upchannelization rate of 4)</summary>

![synthetic_test_2_c4_beam0](./plots/synthetic_test_2_c4_beam0.png)
![synthetic_test_2_c4_beam0_zoom](./plots/synthetic_test_2_c4_beam0_zoom.png)
</details>


</details>