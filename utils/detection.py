import numpy as np
from .cfar_caso_range import cfar_caso_range, CFAR_CASO_Doppler_overlap

def detection(input):
    rangeBinSize = 0.0593
    dopplerFFTSize = 64
    velocityBinSize = 0.0551
    numAntenna = 192
    TDM_MIMO_numTX = 12
    numRxAnt = 16
    sig_integrate = np.sum(np.abs(input)**2, axis=1) + 1
    print("sig_integrate.shape: ", sig_integrate.shape)
    angleFFTSize = 128
    angleBinSkipLeft = 4
    angleBinSkipRight = 4

    N_obj_Rag, Ind_obj_Rag, noise_obj, CFAR_SNR = cfar_caso_range(sig_integrate)
    N_obj = 0
    Ind_obj = []
    detection_results = []

    if N_obj_Rag > 0:
        N_obj, Ind_obj = CFAR_CASO_Doppler_overlap(Ind_obj_Rag, input, sig_integrate)
        
        # Use aggregate noise estimation from the first pass and apply it to objects confirmed in the second pass
        noise_obj_agg = []
        for i_obj in range(N_obj):
            indx1R = Ind_obj[i_obj, 0]
            indx1D = Ind_obj[i_obj, 1]
            ind2R = np.where(Ind_obj_Rag[:, 0] == indx1R)[0]
            ind2D = np.where(Ind_obj_Rag[ind2R, 1] == indx1D)[0]
            noiseInd = ind2R[ind2D]
            noise_obj_agg.append(noise_obj[noiseInd])

        for i_obj in range(N_obj):
            xind = Ind_obj[i_obj, 0] - 1
            detection_result = {}
            detection_result['rangeInd'] = Ind_obj[i_obj, 0] - 1  # range index
            detection_result['range'] = detection_result['rangeInd'] * rangeBinSize  # range estimation

            dopplerInd = Ind_obj[i_obj, 1] - 1  # Doppler index
            detection_result['dopplerInd_org'] = dopplerInd
            detection_result['dopplerInd'] = dopplerInd

            # velocity estimation
            detection_result['doppler'] = (dopplerInd - dopplerFFTSize / 2) * velocityBinSize
            detection_result['doppler_corr'] = detection_result['doppler']
            detection_result['noise_var'] = noise_obj_agg[i_obj]  # noise variance
            detection_result['bin_val'] = np.reshape(input[xind, Ind_obj[i_obj, 1], :], (numAntenna, 1))  # 2d FFT value for the 4 antennas
            detection_result['estSNR'] = np.sum(np.abs(detection_result['bin_val'])**2) / np.sum(detection_result['noise_var'])

            sig_bin = []
            # Only apply max velocity extension if it is enabled and distance is larger than minDisApplyVmaxExtend
            deltaPhi = 2 * np.pi * (dopplerInd - dopplerFFTSize / 2) / (TDM_MIMO_numTX * dopplerFFTSize)
            sig_bin_org = detection_result['bin_val']
            for i_TX in range(TDM_MIMO_numTX):
                RX_ID = np.arange((i_TX - 1) * numRxAnt, i_TX * numRxAnt)
                sig_bin[RX_ID, :] = sig_bin_org[RX_ID, :] * np.exp(-1j * (i_TX - 1) * deltaPhi)
            
            detection_result['bin_val'] = sig_bin
            detection_result['doppler_corr_overlap'] = detection_result['doppler_corr']
            detection_result['doppler_corr_FFT'] = detection_result['doppler_corr']

            detection_results.append(detection_result)

    return detection_results