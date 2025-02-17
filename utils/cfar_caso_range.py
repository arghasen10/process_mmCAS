import numpy as np


def cfar_caso_range(sig):
    cellNum = 8
    gapNum = 8
    K0 = [5, 3]
    
    M_samp, N_pul = sig.shape  
    
    N_obj = 0
    Ind_obj = []
    noise_obj = []
    CFAR_SNR = []

    discardCellLeft = 10
    discardCellRight = 20
    
    for k in range(N_pul):
        sigv = sig[:, k]
        vec = sigv[discardCellLeft:M_samp-discardCellRight]
        vecLeft = vec[:(gapNum + cellNum)]
        vecRight = vec[-(gapNum + cellNum):]
        vec = np.concatenate([vecLeft, vec, vecRight])
        
        for j in range(M_samp - discardCellLeft - discardCellRight):
            cellInd = list(range(j - gapNum - cellNum, j - gapNum)) + list(range(j + gapNum + 1, j + gapNum + cellNum + 1))
            cellInda = list(range(j - gapNum - cellNum, j - gapNum))
            cellIndb = list(range(j + gapNum + 1, j + gapNum + cellNum + 1))

            cellave1a = np.sum(vec[cellInda]) / cellNum
            cellave1b = np.sum(vec[cellIndb]) / cellNum
            cellave1 = min(cellave1a, cellave1b)

            if vec[j + gapNum] > K0[0] * cellave1:
                N_obj += 1
                Ind_obj.append([j + discardCellLeft, k])
                noise_obj.append(cellave1)
                CFAR_SNR.append(vec[j + gapNum] / cellave1)

    for i_obj in range(N_obj):
        ind_range = Ind_obj[i_obj][0]
        ind_Dop = Ind_obj[i_obj][1]
        
        if ind_range <= gapNum:
            cellInd = list(range(ind_range + gapNum + 1, ind_range + gapNum + cellNum)) + \
                      list(range(ind_range + gapNum + 1, ind_range + gapNum + cellNum))
        elif ind_range >= M_samp - gapNum - cellNum:
            cellInd = list(range(ind_range - gapNum - cellNum, ind_range - gapNum - 1)) + \
                      list(range(ind_range - gapNum - cellNum, ind_range - gapNum - 1))
        else:
            cellInd = list(range(ind_range - gapNum - cellNum, ind_range - gapNum - 1)) + \
                      list(range(ind_range + gapNum + 1, ind_range + gapNum + cellNum))

    return N_obj, np.array(Ind_obj), np.array(noise_obj), np.array(CFAR_SNR)


def CFAR_CASO_Doppler_overlap(Ind_obj_Rag, sigCpml, sig_integ):
    maxEnable = 0
    cellNum0 = [8, 4]
    gapNum0 = [8, 0]
    cellNum = cellNum0[1]
    gapNum = gapNum0[1]
    K0 = [5, 3]
    K0 = K0[1]
    numAntenna = 192
    rangeNumBins = sig_integ.shape[0]

    # extract the detected points after range detection
    detected_Rag_Cell = np.unique(Ind_obj_Rag[:, 0])
    sig = sig_integ[detected_Rag_Cell, :]

    M_samp = sig.shape[0]
    N_pul = sig.shape[1]

    # for each point under test, gapNum samples on the two sides are excluded
    # from averaging. Left cellNum/2 and right cellNum/2 samples are used for averaging
    gaptot = gapNum + cellNum

    N_obj = 0
    Ind_obj = []
    noise_obj_an = []
    vec = np.zeros(N_pul + gaptot * 2)

    for k in range(M_samp):
        detected_Rag_Cell_i = detected_Rag_Cell[k]
        ind1 = np.where(Ind_obj_Rag[:, 0] == detected_Rag_Cell_i)
        indR = Ind_obj_Rag[ind1, 1]
        
        # extend the left vector by copying the leftmost the rightmost gaptot samples are not detected
        sigv = sig[k, :]
        vec[:gaptot] = sigv[-gaptot:]
        vec[gaptot:N_pul + gaptot] = sigv
        vec[N_pul + gaptot:] = sigv[:gaptot]

        ind_loc_all = []
        ind_loc_Dop = []
        ind_obj_0 = 0
        noiseEst = np.zeros(N_pul)
        
        for j in range(1 + gaptot, N_pul + gaptot):
            cellInd = np.concatenate([np.arange(j - gaptot, j - gapNum - 1), np.arange(j + gapNum + 1, j + gaptot)])
            noiseEst[j - gaptot] = np.sum(vec[cellInd])

        for j in range(1 + gaptot, N_pul + gaptot):
            j0 = j - gaptot
            cellInd = np.concatenate([np.arange(j - gaptot, j - gapNum - 1), np.arange(j + gapNum + 1, j + gaptot)])
            cellInda = np.arange(j - gaptot, j - gapNum - 1)
            cellIndb = np.arange(j + gapNum + 1, j + gaptot)

            cellave1a = np.sum(vec[cellInda]) / cellNum
            cellave1b = np.sum(vec[cellIndb]) / cellNum
            cellave1 = min(cellave1a, cellave1b)

            maxInCell = np.max(vec[cellInd])
            if maxEnable == 1:
                # detect only if it is the maximum within window
                condition = (vec[j] > K0 * cellave1) and (vec[j] > maxInCell)
            else:
                condition = vec[j] > K0 * cellave1

            if condition:
                # check if this detection overlaps with the Doppler detection
                if np.any(indR == j0):
                    # find overlap, declare a detection
                    ind_win = detected_Rag_Cell_i
                    # range index
                    ind_loc_all.append(ind_win)
                    # Doppler index
                    ind_loc_Dop.append(j0)

        if len(ind_loc_all) > 0:
            ind_obj_0 = np.column_stack([ind_loc_all, ind_loc_Dop])
            if len(Ind_obj) == 0:
                Ind_obj = ind_obj_0
            else:
                # following process is to avoid replicated detection points
                ind_obj_0_sum = ind_loc_all + 10000 * np.array(ind_loc_Dop)
                Ind_obj_sum = Ind_obj[:, 0] + 10000 * Ind_obj[:, 1]
                for ii in range(len(ind_loc_all)):
                    if np.sum(Ind_obj_sum == ind_obj_0_sum[ii]) == 0:
                        Ind_obj = np.vstack([Ind_obj, ind_obj_0[ii, :]])

    N_obj = Ind_obj.shape[0]

    # reset the reference window size to range direction
    cellNum = cellNum0[0]
    gapNum = gapNum0[0]
    gaptot = gapNum + cellNum

    # get the noise variance for each antenna
    N_obj_valid = 0
    Ind_obj_valid = []

    for i_obj in range(N_obj):
        ind_range = Ind_obj[i_obj, 0]
        ind_Dop = Ind_obj[i_obj, 1]

        # skip detected points with signal power less than obj.powerThre
        if ind_Dop < sigCpml.shape[1]:
            if np.min(np.abs(sigCpml[ind_range, ind_Dop, :]) ** 2) < 0:
                continue
        else:
            print(f"Warning: ind_Dop ({ind_Dop}) exceeds the valid range.")
        
        if ind_range <= gaptot:
            # on the left boundary, use the right side samples twice
            cellInd = np.concatenate([np.arange(ind_range + gapNum + 1, ind_range + gaptot), np.arange(ind_range + gapNum + 1, ind_range + gaptot)])
        elif ind_range >= rangeNumBins - gaptot + 1:
            # on the right boundary, use the left side samples twice
            cellInd = np.concatenate([np.arange(ind_range - gaptot, ind_range - gapNum - 1), np.arange(ind_range - gaptot, ind_range - gapNum - 1)])
        else:
            cellInd = np.concatenate([np.arange(ind_range - gaptot, ind_range - gapNum - 1), np.arange(ind_range + gapNum + 1, ind_range + gaptot)])

        N_obj_valid += 1
        noise_obj_an.append(np.mean(np.abs(sigCpml[cellInd, ind_Dop, :]) ** 2, axis=0))

        Ind_obj_valid.append(Ind_obj[i_obj, :])

    N_obj = N_obj_valid
    Ind_obj = np.array(Ind_obj_valid)

    return N_obj, Ind_obj, np.array(noise_obj_an)