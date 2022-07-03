Density matrix calculations for quantum cascade detectors
--------------------------------------------------------------
**Quantum cascade detectors (QCDs)** are mid-infrared sensors, where the transport occurs via successive hops between confined electronic states, even without an applied bias. QCDs advantages include GHz frequency bandwidth and theoretical saturation intensity of the order of hundreds of MW/cm2, unique characteristics for laser-based operations in the long infrared. 

You can find here a python code developed to describe the current transport under illumination in a AlGaAs/GaAs **meta-material enhanced QCD** presented in reference [Bigioli A. et al., 2020](https://aip.scitation.org/doi/10.1063/5.0004591), but it can be adapted to different structures.
The calculation is based on a **density matrix** (DM) description, and it is constructed on the basis of wave functions localized on either side of the extraction barrier and including the coherent tunnelling as a coupling potential. The transport through non-resonant energy subbands is treated with semiclassical rate equations retrieving the Schottky diode analogy in temperature.
It produces the responsivity (and so the photocurrent) as a function of the temperature and applied electric field. It is also possible to extract the population present in every levels in a given condition. The code is composed of:
* **main program**: where the external conditions such as the light intensity, are specified; the parameters of the quantum structure are entered and the photocurrent is calculated and plotted;
* **functions**: where the density matrix and equilibirum carrier distribution calculations are developed. The ouputs are the coherences and the populations at every levels;
* **parameters**: the parameters concerning the band diagram of the chosen structure (material parameters, scattering times, energies, phonon energies...);
* **constants**: physical constants.

Scattering times and the energies related to the specific band diagram have to be calculated separately (their code and values are not provided here).
  
The proposed model efficiently predicts the responsivity of the patch-antenna QCD of ref. [Bigioli A. et al., 2020](https://aip.scitation.org/doi/10.1063/5.0004591), over a broad range of voltages and temperatures, with an easy-to-use approach. As in the tested devices the external radiation is collected by a meta-material composed of a double-metal patch antenna embedding the detector heterostructure, this light coupling enhancement is accounted as a factor to be multiplied to the photocurrent density. The patch-antenna factor is derived in chapter 3 of the reference [Bigioli A., PhD thesis, 2021](http://www.theses.fr/2021UNIP7104), to be downloaded also in [PhD thesis,2021](https://www.researchgate.net/publication/358378549_Uncooled_unipolar_receivers_for_9_m_wavelength_heterodyne_detection). In those references, the chapter 2 is dedicated to the theory underlying the Python code, please refer to this for further theoretical explanations. If you make use of it, please consider citing this document and this code.

The structure of the code has been inspired by the algorithm proposed for quantum cascade lasers in reference [Pegolotti G., PhD thesis, 2014]( http://www.theses.fr/2014PA077136).



