# High-Speed SERDES for PCIe Gen5

This repository contains the MATLAB simulation files used for the SERDES final project.

## Files

### `FFE.m`

System-level validation of the digital RX Feed-Forward Equalizer.

The simulation includes:

* PRBS31 PAM2 data
* Real S-parameter channel model
* 4-bit ADC quantization
* 41-tap digital RX-FFE
* BER comparison before and after FFE
* Eye diagrams and FFE frequency response

### `generate_ENOB.m`

Evaluates the 4-bit ADC performance using the exported ADC output signals.

The script:

* Samples the ADC outputs at the clock edges
* Reconstructs the 4-bit ADC codes
* Calculates SINAD
* Calculates the ADC Effective Number of Bits (ENOB)

### `ADC_FFE_CDR_System_simulation.m`

Top-level ADC and RX-FFE simulation using the actual outputs of four 8-GS/s ADCs.

The four ADC lanes are interleaved to create a 32-GS/s data stream, which is then processed by the 41-tap digital RX-FFE.

The simulation evaluates:

* ADC interleaving
* PRBS31 alignment
* FFE training
* BER before and after FFE
* Digital eye opening
* FFE frequency response

## Project

High-Speed SERDES Receiver for PCIe Gen5-class operation.
