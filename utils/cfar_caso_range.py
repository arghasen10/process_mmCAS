import numpy as np


def cfar_caso_range(sig):
    cellNum = (8,4)
    gapNum = (8,0)
    cellNum = cellNum[0]
    gapNum = gapNum[0]
    K0 = (5,3)

    M_samp=sig.shape[0]
    N_pul=sig.shape[1]
    gaptot=gapNum + cellNum
    N_obj=0
    Ind_obj=[]
    noise_obj = []
    CFAR_SNR = []

    discardCellLeft = 10
    discardCellRight = 20

    for k in range(N_pul):
        sigv = sig[:,k].T
        vec = sigv[discardCellLeft:M_samp-discardCellRight]
        vecLeft = vec[:gaptot]
        vecRight = vec[-gaptot:]
        vec = np.concatenate([vecLeft, vec, vecRight])
        for j in range(M_samp-discardCellLeft-discardCellRight):
            cellInd = list(range(j - gaptot, j - gapNum)) + list(range(j + gapNum + 1, j + gaptot + 1))
            cellInd = np.array(cellInd) + gaptot

            cellInda = np.array(range(j - gaptot, j - gapNum)) + gaptot
            cellIndb = np.array(range(j + gapNum + 1, j + gaptot + 1)) + gaptot

            cellave1a = np.sum(vec[cellInda]) / cellNum
            cellave1b = np.sum(vec[cellIndb]) / cellNum
            cellave1 = min(cellave1a, cellave1b)
            if vec[j + gaptot] > K0 * cellave1:
                N_obj += 1
                Ind_obj.append([j + discardCellLeft, k])
                noise_obj.append(cellave1)  # Save the noise level
                CFAR_SNR.append(vec[j + gaptot] / cellave1)

    for i_obj in range(N_obj):
        ind_range = Ind_obj[i_obj][0]
        ind_Dop = Ind_obj[i_obj][1]

        if ind_range <= gaptot:
            # On the left boundary, use the right-side samples twice
            cellInd = list(range(ind_range + gapNum + 1, ind_range + gaptot + 1)) * 2
        elif ind_range >= M_samp - gaptot + 1:
            # On the right boundary, use the left-side samples twice
            cellInd = list(range(ind_range - gaptot, ind_range - gapNum)) * 2
        else:
            cellInd = (list(range(ind_range - gaptot, ind_range - gapNum)) + 
                    list(range(ind_range + gapNum + 1, ind_range + gaptot + 1)))

    return N_obj, Ind_obj, noise_obj, CFAR_SNR