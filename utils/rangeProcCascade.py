import numpy as np

numAntenna = 16             # number of antennas
numAdcSamplePerChirp = 256  # number of samples per chirp
rangeFFTSize = 256          # FFT size       
rangeWindowEnable  = 1      # flag to enable or disable windowing before range FFT
rangeWindowCoeff  = []      # range FFT window coefficients (one side)
rangeWindowCoeffVec = []    # range FFT window coefficients of length rangeFFTSize
scaleFactorRange  = 1   


def rangeProc(adc_data):
    num_lines = adc_data.shape[1]
    num_ant = adc_data.shape[2]
    out = np.zeros((rangeFFTSize, num_lines, num_ant), dtype=np.complex64)
    
    for i_an in range(num_ant):
        hanning_window = np.hanning(numAdcSamplePerChirp)
        input_mat = np.squeeze(adc_data[:, :, i_an])
        input_mat -= np.mean(input_mat, axis=0)
        input_mat *= hanning_window[:,np.newaxis]
        fft_output = np.fft.fft(input_mat, n=rangeFFTSize, axis=0)
        out[:, :, i_an] = fft_output    
    return out
