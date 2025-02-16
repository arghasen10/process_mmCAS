function main()
    frameIdx = 36;
    frameCountGlobal = frameIdx;
    if ismac
        filepath = '/Users/arghasen/Downloads/cascade_cpu_head_1/';
    elseif ispc
        filepath = 'C:\ti\mmwave_studio_02_01_01_00\mmWaveStudio\PostProc\cascade_cpu_head_1\';
    end
    masterpath = fullfile(filepath, 'master_0000_data.bin');
    slave1path = fullfile(filepath, 'slave1_0000_data.bin');
    slave2path = fullfile(filepath, 'slave2_0000_data.bin');
    slave3path = fullfile(filepath, 'slave3_0000_data.bin');
    [radar_data_Rxchain_master] = readBinFile(masterpath, frameIdx,256,12,64,4,1);
    [radar_data_Rxchain_slave1] = readBinFile(slave1path, frameIdx,256,12,64,4,1);
    [radar_data_Rxchain_slave2] = readBinFile(slave2path, frameIdx,256,12,64,4,1);
    [radar_data_Rxchain_slave3] = readBinFile(slave3path, frameIdx,256,12,64,4,1);
    radar_data_Rxchain(:,:,1:4,:) = radar_data_Rxchain_master;
    radar_data_Rxchain(:,:,5:8,:) = radar_data_Rxchain_slave1;
    radar_data_Rxchain(:,:,9:12,:) = radar_data_Rxchain_slave2;
    radar_data_Rxchain(:,:,13:16,:) = radar_data_Rxchain_slave3;
    rx = [13,14,15,16,1,2,3,4,9,10,11,12,5,6,7,8];
    RxForMIMOProcess = rx;
    radar_data_Rxchain = radar_data_Rxchain(:,:,rx,:); 
    disp(size(radar_data_Rxchain));
    disp("final value" + radar_data_Rxchain(1,1,2,1:10));
    rangeFFTSize = 256;
    rangeFFTOut = [];
    DopplerFFTOut = [];
    cnt = 1;
    [D,IdTxForMIMOProcess] = generateD;
    ind = find(D(:,2)==0);
    [val ID_unique] = unique(D(ind,1));
    antenna_azimuthonly = ind(ID_unique);
    rangeBinSize = 0.0593;
    for i_tx = 1: size(radar_data_Rxchain,4)
        % range FFT
        rangeFFTOut(:,:,:,i_tx)     = rangeProcess(radar_data_Rxchain(:,:,:,i_tx));

        % Doppler FFT
        DopplerFFTOut(:,:,:,i_tx)   = dopplerProcess(rangeFFTOut(:,:,:,i_tx));

    end
    DopplerFFTOut = reshape(DopplerFFTOut,size(DopplerFFTOut,1), size(DopplerFFTOut,2), size(DopplerFFTOut,3)*size(DopplerFFTOut,4));
    sig_integrate = 10*log10(sum((abs(DopplerFFTOut)).^2,3) + 1);
    
    detection_results = detection(DopplerFFTOut);
    detection_results_all{cnt} =  detection_results;

    detect_all_points = [];
    for iobj = 1:length(detection_results)
        detect_all_points (iobj,1)=detection_results(iobj).rangeInd+1;
        detect_all_points (iobj,2)=detection_results(iobj).dopplerInd_org+1;
        detect_all_points (iobj,4)=detection_results(iobj).estSNR;
    end
    
    figure(1);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])                
    subplot(2,2,1)               
    plot((1:size(sig_integrate,1))*rangeBinSize, sig_integrate(:,size(sig_integrate,2)/2+1),'g','LineWidth',4);hold on; grid on
    for ii=1:size(sig_integrate,2)
        plot((1:size(sig_integrate,1))*rangeBinSize, sig_integrate(:,ii));hold on; grid on
        if ~isempty(detection_results)
            ind = find(detect_all_points(:,2)==ii);
            if (~isempty(ind))
                rangeInd = detect_all_points(ind,1);
                plot(rangeInd*rangeBinSize, sig_integrate(rangeInd,ii),'o','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[.49 1 .63],...
                    'MarkerSize',6);
            end
        end
    end

    %title(['FrameID: ' num2str(cnt)]);
    xlabel('Range(m)');
    ylabel('Receive Power (dB)')
    title(['Range Profile(zero Doppler - thick green line): frameID ' num2str(frameIdx)]);
    hold off;
    subplot(2,2,2);
    %subplot_tight(2,2,2,0.1)
    imagesc((sig_integrate))
    c = colorbar;
    c.Label.String = 'Relative Power(dB)';
    title(' Range/Velocity Plot');
    pause(0.01)
    
    angles_all_points = [];
    xyz = [];
    %if 0
    if ~isempty(detection_results)
        % DOA, the results include detection results + angle estimation results.
        % access data with angleEst{frame}(objectIdx).fieldName
        angleEst = DOAEstimate(detection_results);

        numObj = length(detection_results);
        out = detection_results;
        numAoAObjCnt = 0;
        % extended detection_obj to include the angles information
        angleFFTs_detection = zeros(numObj, 256, 256);
        for i_obj = 1:numObj
            current_obj = detection_results(i_obj);
            X = current_obj.bin_val; 
            [angle, angleFFT] = DOA_beamformingFFT_2D(X);
            angleFFTs_detection(i_obj,:,:) = angleFFT;
        end
        angleFFTs_sum = abs(sum(angleFFTs_detection, 1));
        angleFFTs_sum_frames(frameIdx,:,:) = angleFFTs_sum;
        if length(angleEst) > 0
            for iobj = 1:length(angleEst)
                angles_all_points (iobj,1:2)=angleEst(iobj).angles(1:2);
                angles_all_points (iobj,3)=angleEst(iobj).estSNR;
                angles_all_points (iobj,4)=angleEst(iobj).rangeInd;
                angles_all_points (iobj,5)=angleEst(iobj).doppler_corr;
                angles_all_points (iobj,6)=angleEst(iobj).range;
                %switch left and right, the azimuth angle is flipped
                xyz(iobj,1) = angles_all_points (iobj,6)*sind(angles_all_points (iobj,1)*-1)*cosd(angles_all_points (iobj,2));
                xyz(iobj,2) = angles_all_points (iobj,6)*cosd(angles_all_points (iobj,1)*-1)*cosd(angles_all_points (iobj,2));
                %switch upside and down, the elevation angle is flipped
                xyz(iobj,3) = angles_all_points (iobj,6)*sind(angles_all_points (iobj,2)*-1);
                xyz(iobj,4) = angleEst(iobj).doppler_corr;
                xyz(iobj,9) = angleEst(iobj).dopplerInd_org;
                xyz(iobj,5) = angleEst(iobj).range;
                xyz(iobj,6) = angleEst(iobj).estSNR;
                xyz(iobj,7) = angleEst(iobj).doppler_corr_overlap;
                xyz(iobj,8) = angleEst(iobj).doppler_corr_FFT;

            end
            angles_all_all{cnt} = angles_all_points;
            xyz_all{cnt}  = xyz;
            maxRangeShow = rangeBinSize*rangeFFTSize;
            %tic
            moveID = find(abs(xyz(:,4))>=0);
            subplot(2,2,4);                        

            if cnt==1
                scatter3(xyz(moveID,1),xyz(moveID,2),xyz(moveID,3),45,(xyz(moveID,4)),'filled');
            else
                yz = [xyz_all{cnt}; xyz_all{cnt-1}];
                scatter3(xyz(moveID,1),xyz(moveID,2),xyz(moveID,3),45,(xyz(moveID,4)),'filled');
            end

            c = colorbar;
            c.Label.String = 'velocity (m/s)';                        
            grid on;

            xlim([-20 20])
            ylim([1 maxRangeShow])
            %zlim([-4 4])
            zlim([-5 5])
            xlabel('X (m)')
            ylabel('y (m)')
            zlabel('Z (m)')                        

            view([-9 15])                        
            title(' 3D point cloud');

            %plot range and azimuth heatmap
            subplot(2,2,3)
            STATIC_ONLY = 1;
            minRangeBinKeep =  5;
            rightRangeBinDiscard =  20;
            [mag_data_static(:,:,frameCountGlobal) mag_data_dynamic(:,:,frameCountGlobal) y_axis x_axis]= plot_range_azimuth_2D(rangeBinSize, DopplerFFTOut,...
                length(IdTxForMIMOProcess),length(RxForMIMOProcess), ...
                antenna_azimuthonly, 1, STATIC_ONLY, 1, minRangeBinKeep, rightRangeBinDiscard);
            title('range/azimuth heat map static objects')
        end
    end
            
            
    disp("DopplerFFTOut size");
    disp(size(DopplerFFTOut));
    figure;
    d =  abs(sum(DopplerFFTOut,3));
    imagesc(d);
    colormap(jet);
    colorbar;
    disp(d(1,1:10));
    disp('Data Loaded');
end

function [adcData1Complex] = readBinFile(fileFullPath, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices)
    Expected_Num_SamplesPerFrame = numSamplePerChirp*numChirpPerLoop*numLoops*numRXPerDevice*2;
    fp = fopen(fileFullPath, 'r');
    fseek(fp,(frameIdx-1)*Expected_Num_SamplesPerFrame*2, 'bof');
    adcData1 = fread(fp,Expected_Num_SamplesPerFrame,'uint16');
    disp("After file read: " + adcData1(1:10));
    neg = logical(bitget(adcData1, 16));
    disp("neg: " + neg(1:10));
    adcData1(neg) = adcData1(neg) - 2^16;
    disp("After neg operation : " + adcData1(1:10));
    adcData1 = adcData1(1:2:end) + sqrt(-1)*adcData1(2:2:end);
    disp(adcData1(1:11));
    adcData1Complex = reshape(adcData1, numRXPerDevice, numSamplePerChirp, numChirpPerLoop, numLoops);
    adcData1Complex = permute(adcData1Complex, [2 4 1 3]);
    disp(size(adcData1Complex));
    disp("after  " + adcData1Complex(1,1,1,1:11));
    
   %% disp(adcData1Complex(1,1,1,1:11));
    fclose(fp);
end

function [out] = rangeProcess(input)
    windowCoeff = hanning(256);
    rangeWindowCoeff = windowCoeff(1:128);
    rangeWindowCoeffVec         = ones(256, 1);
    rangeWindowCoeffVec(1:128) = rangeWindowCoeff;
    rangeWindowCoeffVec(129:256) = rangeWindowCoeffVec(128:-1:1);
    numLines  = size(input,2);
    numAnt    = size(input,3);
    out = zeros(256, numLines, numAnt);
    for i_an = 1:numAnt  
        inputMat    = squeeze(input(:,:,i_an));
        inputMat    = bsxfun(@minus, inputMat, mean(inputMat));
        inputMat    = bsxfun(@times, inputMat, rangeWindowCoeffVec);
        fftOutput   = fft(inputMat, 256);

        % apply Doppler windowing and scaling to the output. Doppler windowing moved to DopplerProc.m (5/24)
        % fftOutput   = bsxfun(@times, fftOutput, obj.scaleFactorRange*obj.dopplerWindowCoeffVec.');

        % populate in the data cube
        out(:,:,i_an) = fftOutput;

    end            


end


function [out] = dopplerProcess(input)
    numLines  = size(input,1);
    numAnt    = size(input,3);
    numChirpsPerVirAnt = 64;
    dopplerFFTSize = 64;
    out = zeros(numLines, dopplerFFTSize, numAnt);
    windowCoeff = hanning(numChirpsPerVirAnt);
    dopplerWindowCoeff = windowCoeff(1:(numChirpsPerVirAnt/2));
    dopplerWinLen               = length(dopplerWindowCoeff);
    dopplerWindowCoeffVec       = ones(numChirpsPerVirAnt, 1);
    dopplerWindowCoeffVec(1:dopplerWinLen) = dopplerWindowCoeff;
    dopplerWindowCoeffVec(numChirpsPerVirAnt-dopplerWinLen+1:numChirpsPerVirAnt) = dopplerWindowCoeffVec(dopplerWinLen:-1:1);
    for i_an = 1:numAnt                  
                   
        %% vectorized version
        inputMat    = squeeze(input(:,:,i_an));
        inputMat    = bsxfun(@times, inputMat, dopplerWindowCoeffVec.');
        fftOutput   = fft(inputMat, dopplerFFTSize, 2);                         
        fftOutput   = fftshift(fftOutput, 2);
        % populate in the data cube
        out(:,:,i_an) = fftOutput;

    end
    
end

function [detection_results] = detection(input)
    rangeBinSize = 0.0593;
    dopplerFFTSize = 64;
    velocityBinSize = 0.0551;
    numAntenna = 192;
    TDM_MIMO_numTX = 12;
    numRxAnt = 16;
    sig_integrate = sum((abs(input)).^2,3) + 1;
    angleFFTSize = 128;
    angleBinSkipLeft = 4;
    angleBinSkipRight = 4;
    [N_obj_Rag, Ind_obj_Rag, noise_obj, CFAR_SNR] = CFAR_CASO_Range(sig_integrate);
    N_obj = 0;
    Ind_obj = [];
    detection_results = {};
    if (N_obj_Rag>0)
        [N_obj, Ind_obj] = CFAR_CASO_Doppler_overlap(Ind_obj_Rag, input, sig_integrate);
        detection_results = [];
        
        % Use aggregate noise estimation from the first pass and apply
        % it to objects confirmed in the second pass
        noise_obj_agg = [];
        for i_obj = 1:N_obj
            indx1R = Ind_obj(i_obj,1);
            indx1D = Ind_obj(i_obj,2);
            ind2R = find(Ind_obj_Rag(:,1) == indx1R);
            ind2D = find(Ind_obj_Rag(ind2R,2) == indx1D);
            noiseInd = ind2R(ind2D);
            noise_obj_agg(i_obj) = noise_obj(noiseInd);
        end
        
        for i_obj = 1:N_obj
            xind = (Ind_obj(i_obj,1)-1) +1;
            detection_results(i_obj).rangeInd = Ind_obj(i_obj, 1) - 1;  %range index
            detection_results(i_obj).range = (detection_results(i_obj).rangeInd) * rangeBinSize;  %range estimation
            dopplerInd  = Ind_obj(i_obj, 2) - 1;  %Doppler index
            detection_results(i_obj).dopplerInd_org = dopplerInd;
            detection_results(i_obj).dopplerInd = dopplerInd;
            
            %velocity estimation
            detection_results(i_obj).doppler = (dopplerInd-dopplerFFTSize/2)*velocityBinSize;
            detection_results(i_obj).doppler_corr = detection_results (i_obj).doppler;
            detection_results(i_obj).noise_var = noise_obj_agg(i_obj);       %noise variance
            detection_results(i_obj).bin_val  = reshape(input(xind, Ind_obj(i_obj,2),:),numAntenna,1);  %2d FFT value for the 4 antennas
            %detection_results(i_obj).estSNR  = 10*log10(sum(abs(detection_results (i_obj).bin_val).^2)/sum(detection_results (i_obj).noise_var));  %2d FFT value for the 4 antennas
            detection_results(i_obj).estSNR  = (sum(abs(detection_results (i_obj).bin_val).^2)/sum(detection_results (i_obj).noise_var));  
            
            sig_bin = [];
            %only apply max velocity extention if it is enabled and distance is larger
            %than minDisApplyVmaxExtend
            
            deltaPhi = 2*pi*(dopplerInd-dopplerFFTSize/2)/( TDM_MIMO_numTX*dopplerFFTSize);
            sig_bin_org = detection_results (i_obj).bin_val;
            for i_TX = 1:TDM_MIMO_numTX
                RX_ID = (i_TX-1)*numRxAnt+1 : i_TX*numRxAnt;
                sig_bin(RX_ID,: )= sig_bin_org(RX_ID )* exp(-1j*(i_TX-1)*deltaPhi);
            end
            detection_results(i_obj).bin_val = sig_bin;
            detection_results(i_obj).doppler_corr_overlap = detection_results(i_obj).doppler_corr;
            detection_results(i_obj).doppler_corr_FFT = detection_results(i_obj).doppler_corr;
        end
    end
end

function [N_obj, Ind_obj, noise_obj, CFAR_SNR] = CFAR_CASO_Range(sig)

    cellNum = [8,4];
    gapNum = [8,0];
    cellNum = cellNum(1);
    gapNum = gapNum(1);
    K0 = [5,3];

    M_samp=size(sig, 1);
    N_pul=size(sig, 2);
    %% input:  assuming size(input) = [numSamplePerChipr numChirpsPerFrame numAntenna]


    %for each point under test, gapNum samples on the two sides are excluded
    %from averaging. Left cellNum/2 and right cellNum/2 samples are used for
    %averaging
    gaptot=gapNum + cellNum;
    N_obj=0;
    Ind_obj=[];
    noise_obj = [];
    CFAR_SNR = [];

    discardCellLeft = 10;
    discardCellRight = 20;


    %for the first gaptot samples only use the right sample
    for k=1:N_pul
        sigv=(sig(:,k))';
        vec = sigv(discardCellLeft+1:M_samp-discardCellRight);
        vecLeft = vec(1:(gaptot));
        vecRight = vec(end-(gaptot)+1:end);
        vec = [vecLeft vec vecRight];
        for j=1:(M_samp-discardCellLeft-discardCellRight)
            cellInd=[j-gaptot: j-gapNum-1 j+gapNum+1:j+gaptot];
            cellInd=cellInd + gaptot;
            cellInda=[j-gaptot: j-gapNum-1];
            cellInda=cellInda + gaptot;
            cellIndb=[ j+gapNum+1:j+gaptot];
            cellIndb=cellIndb + gaptot;

            cellave1a =sum(vec(cellInda))/(cellNum);
            cellave1b =sum(vec(cellIndb))/(cellNum);
            cellave1 = min(cellave1a,cellave1b);
            if vec(j+gaptot)>K0*cellave1
                N_obj=N_obj+1;
                Ind_obj(N_obj,:)=[j+discardCellLeft, k];
                noise_obj(N_obj) = cellave1; %save the noise level
                CFAR_SNR(N_obj) = vec(j+gaptot)/cellave1;
            end        
        end
    end

    %get the noise variance for each antenna
    for i_obj = 1:N_obj

        ind_range = Ind_obj(i_obj,1);
        ind_Dop = Ind_obj(i_obj,2);
        if ind_range<= gaptot
            %on the left boundary, use the right side samples twice
            cellInd=[ind_range+gapNum+1:ind_range+gaptot ind_range+gapNum+1:ind_range+gaptot];
        elseif ind_range>=M_samp-gaptot+1
            %on the right boundary, use the left side samples twice
            cellInd=[ind_range-gaptot: ind_range-gapNum-1 ind_range-gaptot: ind_range-gapNum-1];
        else
            cellInd=[ind_range-gaptot: ind_range-gapNum-1 ind_range+gapNum+1:ind_range+gaptot];

        end


    end
end




function [N_obj, Ind_obj, noise_obj_an] = CFAR_CASO_Doppler_overlap(Ind_obj_Rag, sigCpml, sig_integ)
    maxEnable = 0;
    cellNum0 = [8,4];
    gapNum0 = [8,0];
    cellNum = cellNum0(2);
    gapNum = gapNum0(2);
    K0 = [5, 3];
    K0 = K0(2);
    numAntenna = 192;
    rangeNumBins = size(sig_integ,1);

    %extract the detected points after range detection
    detected_Rag_Cell = unique(Ind_obj_Rag(:,1));
    sig = sig_integ(detected_Rag_Cell,:);

    M_samp=size(sig, 1);
    N_pul=size(sig, 2);


    %for each point under test, gapNum samples on the two sides are excluded
    %from averaging. Left cellNum/2 and right cellNum/2 samples are used for
    %averaging
    gaptot=gapNum + cellNum;

    N_obj=0;
    Ind_obj=[];
    noise_obj_an = [];
    vec=zeros(1,N_pul+gaptot*2);
    for k=1:M_samp
        %get the range index at current range index
        detected_Rag_Cell_i = detected_Rag_Cell(k);
        ind1 = find(Ind_obj_Rag(:,1) == detected_Rag_Cell_i);
        indR = Ind_obj_Rag(ind1, 2);
        %extend the left the vector by copying the left most the right most
        %gaptot samples are not detected.
        sigv=(sig(k,:));
        vec(1:gaptot) = sigv(end-gaptot+1:end);
        vec(gaptot+1: N_pul+gaptot) = sigv;
        vec(N_pul+gaptot+1:end) = sigv(1:gaptot);
        %start to process
        ind_loc_all = [];
        ind_loc_Dop = [];
        ind_obj_0 = 0;
        noiseEst = zeros(1,N_pul);
        for j=1+gaptot:N_pul+gaptot
            cellInd=[j-gaptot: j-gapNum-1 j+gapNum+1:j+gaptot];
            noiseEst(j-gaptot) = sum(vec(cellInd));
        end
        for j=1+gaptot:N_pul+gaptot
            j0 = j - gaptot;
            cellInd=[j-gaptot: j-gapNum-1 j+gapNum+1:j+gaptot];
            cellInda = [j-gaptot: j-gapNum-1 ];
            cellIndb =[j+gapNum+1:j+gaptot];

            cellave1a =sum(vec(cellInda))/(cellNum);
            cellave1b =sum(vec(cellIndb))/(cellNum);
            cellave1 = min(cellave1a,cellave1b);        

            maxInCell = max(vec(cellInd));
            if maxEnable==1
                %detect only if it is the maximum within window
                condition = ((vec(j)>K0*cellave1)) && ((vec(j)>maxInCell));
            else
                condition = vec(j)>K0*cellave1;
            end

            if condition==1
                %check if this detection overlap with the Doppler detection
                if(find(indR == j0))
                    %find overlap, declare a detection
                    ind_win = detected_Rag_Cell_i;
                    %range index
                    ind_loc_all = [ind_loc_all ind_win];
                    %Doppler index
                    ind_loc_Dop = [ind_loc_Dop j0];
                end

            end

        end
        ind_obj_0 = [];


        if (length(ind_loc_all)>0)
            ind_obj_0(:,1) = ((ind_loc_all));
            ind_obj_0(:,2) = ind_loc_Dop;
            if size(Ind_obj,1) ==0
                Ind_obj = ind_obj_0;
            else

                %following process is to avoid replicated detection points
                ind_obj_0_sum = ind_loc_all + 10000*ind_loc_Dop;
                Ind_obj_sum = Ind_obj(:,1) + 10000*Ind_obj(:,2);
                for ii= 1: length(ind_loc_all)
                    if (length(find(Ind_obj_sum == ind_obj_0_sum(ii)))==0)
                        Ind_obj = [Ind_obj ; ind_obj_0(ii,:)];
                    end
                end
            end
        end

    end

    N_obj = size(Ind_obj,1);

    %reset the ref window size to range direction
    cellNum = cellNum0(1);
    gapNum = gapNum0(1);
    gaptot=gapNum + cellNum;
    %get the noise variance for each antenna
    N_obj_valid = 0;
    Ind_obj_valid = [];
    for i_obj = 1:N_obj    
        ind_range = Ind_obj(i_obj,1);
        ind_Dop = Ind_obj(i_obj,2);
        %skip detected points with signal power less than obj.powerThre
        if (min(abs(sigCpml(ind_range, ind_Dop,:)).^2) < 0)
            continue;
        end
        if ind_range<= gaptot
            %on the left boundary, use the right side samples twice
            cellInd=[ind_range+gapNum+1:ind_range+gaptot ind_range+gapNum+1:ind_range+gaptot];
        elseif ind_range>=rangeNumBins-gaptot+1
            %on the right boundary, use the left side samples twice
            cellInd=[ind_range-gaptot: ind_range-gapNum-1 ind_range-gaptot: ind_range-gapNum-1];
        else
            cellInd=[ind_range-gaptot: ind_range-gapNum-1 ind_range+gapNum+1:ind_range+gaptot];

        end

        N_obj_valid = N_obj_valid +1;
        noise_obj_an(:, i_obj) = reshape((mean(abs(sigCpml(cellInd, ind_Dop, :)).^2, 1)), numAntenna, 1, 1);

        Ind_obj_valid(N_obj_valid,:) = Ind_obj(i_obj,:);    

    end

    N_obj = N_obj_valid;
    Ind_obj = Ind_obj_valid;
end


function [out] = DOAEstimate(detection_results)
    numObj = length(detection_results);
    out = detection_results;
    numAoAObjCnt = 0;
    % extended detection_obj to include the angles information
    for i_obj = 1:numObj
        current_obj = detection_results(i_obj);
        estSNR = 10*log10(sum(abs(current_obj.bin_val).^2)/sum(current_obj.noise_var));
        X = current_obj.bin_val; 
        R = X*X';
        %2D beamforming angle estimation, azimuth is estimated based on 1D FFT output                        
        [DOA_angles angle_sepc_2D_fft ]= DOA_beamformingFFT_2D(X);
        if (numAoAObjCnt == 0)
            out = [];
        end

        for i_obj = 1:size(DOA_angles,2)
            numAoAObjCnt = numAoAObjCnt+1;
            out(numAoAObjCnt).rangeInd = current_obj.rangeInd;
            out(numAoAObjCnt).dopplerInd = current_obj.dopplerInd;
            out(numAoAObjCnt).range = current_obj.range;
            out(numAoAObjCnt).doppler_corr = current_obj.doppler_corr;
            out(numAoAObjCnt).dopplerInd_org = current_obj.dopplerInd_org;

            out(numAoAObjCnt).noise_var = current_obj.noise_var;
            out(numAoAObjCnt).bin_val = current_obj.bin_val;
            out(numAoAObjCnt).estSNR = current_obj.estSNR;
            out(numAoAObjCnt).doppler_corr_overlap = current_obj.doppler_corr_overlap;
            out(numAoAObjCnt).doppler_corr_FFT = current_obj.doppler_corr_FFT;

            out(numAoAObjCnt).angles = DOA_angles(:,i_obj);
            out(numAoAObjCnt).spectrum = angle_sepc_2D_fft;


        end
    end
end


function [angleObj_est angle_sepc_2D_fft]= DOA_beamformingFFT_2D(sig)
    %field of view to do beamforming
    angles_DOA_az = [-80, 80];
    angles_DOA_ele = [-20,20];

    %distance unit in terms of wavelength
    d = 0.5095;
    %2D matrix providing antenna coordinates
    [D,IdTxForMIMOProcess] = generateD;
    DOAFFTSize = 256;
    angleFFTSize = DOAFFTSize;


    %FFT based implementation
    %first form a 2D matrix based on the antenna coordinates
    D = D + 1;
    apertureLen_azim = max(D(:,1));
    apertureLen_elev = max(D(:,2));
    sig_2D = zeros(apertureLen_azim,apertureLen_elev);
    for i_line = 1:apertureLen_elev
        ind = find(D(:,2) == i_line);
        D_sel = D(ind,1);
        sig_sel = sig(ind);
        [val indU] = unique(D_sel);

        sig_2D(D_sel(indU),i_line) = sig_sel(indU);

    end

    %run FFT on azimuth and elevation
    angle_sepc_1D_fft=fftshift(fft(sig_2D,angleFFTSize,1),1); 
    angle_sepc_2D_fft=fftshift(fft(angle_sepc_1D_fft,angleFFTSize,2),2); 

    wx_vec=[-pi:2*pi/angleFFTSize:pi];
    wz_vec=[-pi:2*pi/angleFFTSize:pi];
    wx_vec = wx_vec(1:end-1);
    wz_vec = wz_vec(1:end-1);
    %use one row with complete azimuth antenna of 1D FFT output for azimuth
    %estimation
    spec_azim = abs(angle_sepc_1D_fft(:,1));
    sidelobeLevel_dB_azim = 1;
    sidelobeLevel_dB = sidelobeLevel_dB_azim;
    [peakVal_azim, peakLoc_azim] = DOA_BF_PeakDet_loc(spec_azim);

    if apertureLen_elev ==1
        %azimuth array only, no elevation antennas
        obj_cnt = 1;
        angleObj_est= [];
        for i_obj = 1:length(peakLoc_azim)
            ind = peakLoc_azim(i_obj);

            azim_est = asind(wx_vec(ind)/(2*pi*d));
            if (azim_est >= angles_DOA_az(1) && azim_est <= angles_DOA_az(2))
                angleObj_est(1,obj_cnt) = azim_est;
                angleObj_est(2,obj_cnt) = 0;
                angleObj_est(3,obj_cnt) = ind;
                angleObj_est(4,obj_cnt) = 0;
                obj_cnt = obj_cnt+1;

            else
                continue;
            end
        end

    else
        %azimuth and elevation angle estimation

        % figure(1);plot(spec_azim); hold on; grid on
        % plot(peakLoc_azim, spec_azim(peakLoc_azim),'ro');hold on

        %for each detected azimuth, estimate its elevation
        % figure(2)
        obj_cnt = 1;
        angleObj_est= [];
        sidelobeLevel_dB_elev = 0;
        sidelobeLevel_dB = sidelobeLevel_dB_elev;
        for i_obj = 1:length(peakLoc_azim)
            ind = peakLoc_azim(i_obj);
            spec_elev = abs(angle_sepc_2D_fft(ind,:));
            [peakVal_elev, peakLoc_elev] = DOA_BF_PeakDet_loc(spec_elev);
            %calcualte the angle values
            for j_elev = 1:length(peakVal_elev)
                azim_est = asind(wx_vec(ind)/(2*pi*d));
                elev_est = asind(wz_vec(peakLoc_elev(j_elev))/(2*pi*d));

                if (azim_est >= angles_DOA_az(1) && azim_est <= angles_DOA_az(2) ...
                        &&elev_est >= angles_DOA_ele(1) && elev_est <= angles_DOA_ele(2))
                    angleObj_est(1,obj_cnt) = azim_est;
                    angleObj_est(2,obj_cnt) = elev_est;              
                    angleObj_est(3,obj_cnt) = ind;
                    angleObj_est(4,obj_cnt) = peakLoc_elev(j_elev);
                    %plot(angleObj_est(4,obj_cnt),angleObj_est(3,obj_cnt) ,'x','MarkerSize',12, 'LineWidth',2);
                    %hold on
                    obj_cnt = obj_cnt+1;

                else
                    continue;
                end
            end        
        end    
        %hold off

    end
end

function [D, IdTxForMIMOProcess] = generateD()
    TxToEnable = [12  11  10   9   8   7   6   5   4   3   2   1];
    TI_Cascade_RX_ID = [13 14 15 16 1 2 3 4 9 10 11 12 5 6 7 8 ];
    TI_Cascade_RX_position_azi = [ 11:14 50:53 46:49 0:3  ];
    TxForMIMOProcess = TxToEnable;
    TI_Cascade_TX_position_azi = [11 10 9 32 28 24 20 16 12 8 4 0 ];%12 TX antenna azimuth position on TI 4-chip cascade EVM
    TI_Cascade_TX_position_ele = [6 4 1 0 0 0 0 0 0 0 0 0];%12 TX antenna elevation position on TI 4-chip cascade EVM
    TI_Cascade_RX_position_ele = zeros(1,16);%16 RX antenna elevation position on TI 4-chip cascade EVM
    [IdTxForMIMOProcess ia ib] = intersect(TxForMIMOProcess, TxToEnable,'stable' );
    RxForMIMOProcess = TI_Cascade_RX_ID; %using all 16 RXs, user can also choose subset of RXs for MIMO data analysis
    D_TX = TI_Cascade_TX_position_azi(TxToEnable(ib)); %TX azimuth antenna coordinates
    D_TX_ele = TI_Cascade_TX_position_ele(TxToEnable(ib));%TX elevation antenna coordinates
    D_RX = TI_Cascade_RX_position_azi(RxForMIMOProcess); %RX azimuth antenna coordinate
    D_RX_ele = TI_Cascade_RX_position_ele(RxForMIMOProcess);%RX elevation antenna coordinate
    plotArray = 0;
    RX_id_tot = [];
    RX_id_tot_ele = [];
    for ii = 1:length(D_TX)
        RX_id_new = D_RX + sum(D_TX(ii));
        RX_id_tot = [RX_id_tot RX_id_new];
        RX_id_new_ele = D_RX_ele + D_TX_ele(ii);
        RX_id_tot_ele = [RX_id_tot_ele RX_id_new_ele];
    end
    D(:,1) = RX_id_tot;
    D(:,2) = RX_id_tot_ele;
end

function [peakVal, peakLoc] = DOA_BF_PeakDet_loc(inData)
    gamma = 1.0471;
    sidelobeLevel_dB = 0;

    inData = inData(:);

    minVal = Inf;
    maxVal = 0;
    maxLoc = 0;
    maxData = [];

    locateMax = 0;  % at beginning, not ready for peak detection
    km = 1;         % constant value used in variance calculation

    numMax = 0;
    extendLoc = 0;
    initStage = 1;
    absMaxValue = 0;

    i = 0;
    N = length(inData);
    while (i < (N + extendLoc - 1))
        i = i+1;
        i_loc = rem(i-1, N)+1;
        currentVal = inData(i_loc);
        % record the maximum value
        if currentVal > absMaxValue
            absMaxValue = currentVal;
        end
        % record the current max value and location
        if currentVal > maxVal
            maxVal = currentVal;
            maxLoc = i_loc;
            maxLoc_r = i;
        end

        % record for the current min value and location
        if currentVal < minVal,
            minVal = currentVal;
        end

        if locateMax
            if currentVal < maxVal/gamma
                numMax = numMax + 1;
                bwidth = i - maxLoc_r;
                % Assign maximum value only if the value has fallen below the max by
                % gamma, thereby declaring that the max value was a peak
                maxData = [maxData(1:numMax-1,:) ; maxLoc maxVal bwidth maxLoc_r];

                minVal = currentVal;
                locateMax = 0;
            end
        else
            if currentVal > minVal*gamma
                % Assign minimum value if the value has risen above the min by
                % gamma, thereby declaring that the min value was a valley
                locateMax = 1;
                maxVal = currentVal;
                if (initStage == 1)
                    extendLoc = i;
                    initStage = 0;
                end
            end
        end
    end


    %make sure the max value needs to be cetain dB higher than the side lobes
    %to declare any detection
    estVar = zeros(numMax, 1);
    peakVal = zeros(numMax, 1);
    peakLoc = zeros(numMax, 1);
    delta = [];

    %[v ind] = max(maxData(:,2));
    %peakMean = mean(maxData([1:(ind-1) ind+1:end],2));
    %SNR_DOA = v/peakMean;
    %if v>peakMean*maxPeakThre
    % if the max is different by more than sidelobeLevel_dB dB from the
    % peak, then removed it as a sidelobe
    maxData_ = [];
    numMax_ = 0;
    totPower = 0;
    for i = 1:numMax
        if maxData(i, 2) >= absMaxValue * (10^(-sidelobeLevel_dB/10))
            numMax_ = numMax_ + 1;
            maxData_(numMax_,:) = maxData(i, :);
            totPower = totPower + maxData(i, 2);
        end
    end
    maxData = maxData_;
    numMax = numMax_;

    estVar = zeros(numMax, 1);
    peakVal = zeros(numMax, 1);
    peakLoc = zeros(numMax, 1);

    delta = [];
    for ind = 1:numMax
        peakVal(ind) = maxData(ind,2);
        peakLoc(ind) = rem(maxData(ind,1)-1, N)+1;
    end
end


function  [mag_data_static mag_data_dynamic y_axis x_axis] = plot_range_azimuth_2D(range_resolution, radar_data_pre_3dfft,TDM_MIMO_numTX,numRxAnt,antenna_azimuthonly, LOG, STATIC_ONLY, PLOT_ON, minRangeBinKeep,  rightRangeBinDiscard)
    dopplerFFTSize = size(radar_data_pre_3dfft,2);
    rangeFFTSize = size(radar_data_pre_3dfft,1);
    angleFFTSize = 256;
    % ratio used to decide engergy threshold used to pick non-zero Doppler bins
    ratio = 0.5;
    DopplerCorrection = 0;
    
    if DopplerCorrection == 1
        % add Doppler correction before generating the heatmap
        radar_data_pre_3dfft_DopCor= [];
        for dopplerInd = 1: dopplerFFTSize
            deltaPhi = 2*pi*(dopplerInd-1-dopplerFFTSize/2)/( TDM_MIMO_numTX*dopplerFFTSize);
            sig_bin_org =squeeze(radar_data_pre_3dfft(:,dopplerInd,:));
            for i_TX = 1:TDM_MIMO_numTX
                RX_ID = (i_TX-1)*numRxAnt+1 : i_TX*numRxAnt;
                corVec = repmat(exp(-1j*(i_TX-1)*deltaPhi), rangeFFTSize, numRxAnt);
                radar_data_pre_3dfft_DopCor(:,dopplerInd, RX_ID)= sig_bin_org(:,RX_ID ).* corVec;
            end
        end
        
        radar_data_pre_3dfft = radar_data_pre_3dfft_DopCor;
    end
    radar_data_pre_3dfft = radar_data_pre_3dfft(:,:,antenna_azimuthonly);
    
    radar_data_angle_range = fft(radar_data_pre_3dfft, angleFFTSize, 3);
    n_angle_fft_size = size(radar_data_angle_range,3);
    n_range_fft_size = size(radar_data_angle_range,1);
    
    
    %decide non-zerp doppler bins to be used for dynamic range-azimuth heatmap
    DopplerPower = sum(mean((abs(radar_data_pre_3dfft(:,:,:))),3),1);
    DopplerPower_noDC = DopplerPower([1: dopplerFFTSize/2-1 dopplerFFTSize/2+3:end]);
    [peakVal peakInd] = max(DopplerPower_noDC);
    threshold = peakVal*ratio;
    indSel = find(DopplerPower_noDC >threshold);
    for ii = 1:length(indSel)
        if indSel(ii) > dopplerFFTSize/2-1
            indSel(ii) = indSel(ii) + 3;
        end
    end
    
    radar_data_angle_range_dynamic = squeeze(sum(abs(radar_data_angle_range(:,indSel,:)),2));
    radar_data_angle_range_Static = squeeze(sum(abs(radar_data_angle_range(:,dopplerFFTSize/2+1,:)),2));
    
    
    indices_1D = (minRangeBinKeep:n_range_fft_size-rightRangeBinDiscard);
    max_range = (n_range_fft_size-1)*range_resolution;
    max_range = max_range/2;
    d = 1;
    
    %generate range/angleFFT for zeroDoppler and non-zero Doppler respectively
    radar_data_angle_range_dynamic = fftshift(radar_data_angle_range_dynamic,2);
    radar_data_angle_range_Static = fftshift(radar_data_angle_range_Static,2);
    
    
    sine_theta = -2*((-n_angle_fft_size/2:n_angle_fft_size/2)/n_angle_fft_size)/d;
    cos_theta = sqrt(1-sine_theta.^2);
    
    [R_mat, sine_theta_mat] = meshgrid(indices_1D*range_resolution,sine_theta);
    [~, cos_theta_mat] = meshgrid(indices_1D,cos_theta);
    
    x_axis = R_mat.*cos_theta_mat;
    y_axis = R_mat.*sine_theta_mat;
    mag_data_dynamic = squeeze(abs(radar_data_angle_range_dynamic(indices_1D+1,[1:end 1])));
    mag_data_static = squeeze(abs(radar_data_angle_range_Static(indices_1D+1,[1:end 1])));
    
    
    mag_data_dynamic = mag_data_dynamic';
    mag_data_static = mag_data_static';
    mag_data_dynamic = flipud(mag_data_dynamic);
    mag_data_static = flipud(mag_data_static);
    
    
    if PLOT_ON
        log_plot = LOG;
        if STATIC_ONLY == 1
            if log_plot
                surf(y_axis, x_axis, (mag_data_static).^0.4,'EdgeColor','none');
            else
                surf(y_axis, x_axis, abs(mag_data_static),'EdgeColor','none');
            end
        else
            if log_plot
                surf(y_axis, x_axis, (mag_data_dynamic).^0.4,'EdgeColor','none');
            else
                surf(y_axis, x_axis, abs(mag_data_dynamic),'EdgeColor','none');
            end
        end
        
        view(2);
        xlabel('meters')
        ylabel('meters')
        
    end
end