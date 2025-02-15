function main()
    [radar_data_Rxchain_master] = readBinFile("C:\ti\mmwave_studio_02_01_01_00\mmWaveStudio\PostProc\cascade_cpu_head_1\master_0000_data.bin", 2,256,12,64,4,1);
    [radar_data_Rxchain_slave1] = readBinFile("C:\ti\mmwave_studio_02_01_01_00\mmWaveStudio\PostProc\cascade_cpu_head_1\slave1_0000_data.bin", 2,256,12,64,4,1);
    [radar_data_Rxchain_slave2] = readBinFile("C:\ti\mmwave_studio_02_01_01_00\mmWaveStudio\PostProc\cascade_cpu_head_1\slave2_0000_data.bin", 2,256,12,64,4,1);
    [radar_data_Rxchain_slave3] = readBinFile("C:\ti\mmwave_studio_02_01_01_00\mmWaveStudio\PostProc\cascade_cpu_head_1\slave3_0000_data.bin", 2,256,12,64,4,1);
    radar_data_Rxchain(:,:,1:4,:) = radar_data_Rxchain_master;
    radar_data_Rxchain(:,:,5:8,:) = radar_data_Rxchain_slave1;
    radar_data_Rxchain(:,:,9:12,:) = radar_data_Rxchain_slave2;
    radar_data_Rxchain(:,:,13:16,:) = radar_data_Rxchain_slave3;
    rx = [13    14    15    16     1     2     3     4     9    10    11    12     5     6     7     8];
    radar_data_Rxchain = radar_data_Rxchain(:,:,rx,:); 
    disp(size(radar_data_Rxchain));
    disp("final value" + radar_data_Rxchain(1,1,2,1:10));
    rangeFFTOut = [];
    DopplerFFTOut = [];

    for i_tx = 1: size(radar_data_Rxchain,4)
        % range FFT
        rangeFFTOut(:,:,:,i_tx)     = rangeProcess(radar_data_Rxchain(:,:,:,i_tx));

        % Doppler FFT
        DopplerFFTOut(:,:,:,i_tx)   = datapath(DopplerFFTObj, rangeFFTOut(:,:,:,i_tx));

    end
    figure;
    d =  abs(sum(sum(rangeFFTOut,4),3));
    imagesc(d);
    colormap(jet);
    colorbar;
    disp(d(1:10));
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
    
    out = zeros(numLines, 64, numAnt);
    windowCoeff = hanning(256);
    rangeWindowCoeff = windowCoeff(1:128);
    rangeWindowCoeffVec         = ones(256, 1);
    rangeWindowCoeffVec(1:128) = rangeWindowCoeff;
    rangeWindowCoeffVec(129:256) = rangeWindowCoeffVec(128:-1:1);
    out = zeros(256, numLines, numAnt);
    
    numAntenna = 0           % number of antennas
    numDopplerLines = 0        % number of Doppler lines
    dopplerFFTSize  = 0        % Doppler FFT size
    numChirpsPerVirAnt = 0       
    dopplerWindowEnable  = 0   % flag to enable or disable windowing before doppler FFT
    dopplerWindowCoeff  = []    % Doppler FFT window coefficients (one side)
    dopplerWindowCoeffVec = []  % Doppler FFT window coefficients of length dopplerFFTSize
    scaleFactorDoppler = 0       
    clutterRemove = 0
    FFTOutScaleOn = 0
    
    windowCoeff = hanning(numChirpsPerVirAnt);
    DopplerProcClutterRemove_dopplerWindowCoeff = windowCoeff(1:(64/2));
    DopplerProcClutterRemove_scaleFactorDoppler  = scaleFactor(log2(DopplerFFTSize) - 3);
    DopplerProcClutterRemove_FFTOutScaleOn       = 0; %1: apply scaleFactorRange; 0: scaling factor not applied
    DopplerProcClutterRemove_clutterRemove       = 0;

    dopplerWinLen               = length(DopplerProcClutterRemove_dopplerWindowCoeff);
    dopplerWindowCoeffVec       = ones(64, 1);
    dopplerWindowCoeffVec(1:dopplerWinLen) = obj.dopplerWindowCoeff;
    dopplerWindowCoeffVec(obj.numChirpsPerVirAnt-dopplerWinLen+1:obj.numChirpsPerVirAnt) = dopplerWindowCoeffVec(dopplerWinLen:-1:1);
    obj.dopplerWindowCoeffVec   = dopplerWindowCoeffVec;  
    
    for i_an = 1:numAnt                  
                   
        %% vectorized version
        inputMat    = squeeze(input(:,:,i_an));
        inputMat    = bsxfun(@times, inputMat, obj.dopplerWindowCoeffVec.');
        if obj.clutterRemove ==1
            inputMat = inputMat - (repmat(mean(inputMat'),size(inputMat,2),1))';
        end
        fftOutput   = fft(inputMat, obj.dopplerFFTSize, 2);                         

        if obj.FFTOutScaleOn ==1
            fftOutput   = fftshift(fftOutput, 2) * obj.scaleFactorDoppler;
        else
            fftOutput   = fftshift(fftOutput, 2);
        end
        % populate in the data cube
        out(:,:,i_an) = fftOutput;

    end
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


 