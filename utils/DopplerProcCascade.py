import numpy as np

numAntenna = 16             # number of antennas
numAdcSamplePerChirp = 256  # number of samples per chirp
dopplerFFTSize = 64          # FFT size       
rangeWindowEnable  = 1      # flag to enable or disable windowing before range FFT
rangeWindowCoeff  = []      # range FFT window coefficients (one side)
rangeWindowCoeffVec = []    # range FFT window coefficients of length rangeFFTSize
scaleFactorRange  = 1   


def dopplerProc(adc_data):
    num_lines = adc_data.shape[0]
    num_ant = adc_data.shape[2]
    out = np.zeros((num_lines, dopplerFFTSize, num_ant), dtype=np.complex64)
    
    for i_an in range(num_ant):
        input_mat = np.squeeze(adc_data[:, :, i_an])
        fft_output = np.fft.fft(input_mat, n=dopplerFFTSize, axis=1)
        np.fft.fftshift(fft_output,1)
        out[:, :, i_an] = fft_output    
    return out
