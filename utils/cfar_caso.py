import numpy as np
from .cfar_caso_range import cfar_caso_range

def cfar_case(dopplerFFtOut):
    sig_integrate = np.sum(np.abs(dopplerFFtOut)**2, axis=2) + 1
    angleFFTSize = 128
    angleBinSkipLeft = 4
    angleBinSkipRight = 4
    N_obj, Ind_obj, noise_obj, CFAR_SNR = cfar_caso_range(sig_integrate)
    print(f"N_obj: {N_obj}, Ind_obj: {Ind_obj}, noise_obj: {noise_obj}, CFAR_SNR: {CFAR_SNR}")