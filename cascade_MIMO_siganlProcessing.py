import glob
import numpy as np
from utils.read_ADC_bin_TDA2_separateFiles import read_adc_bin_tda2_separate_files
from utils.rangeProcCascade import rangeProc
from utils.DopplerProcCascade import dopplerProc
from utils.cfar_caso import cfar_case
import seaborn as sns

num_sample_per_chirp=256
num_chirp_per_loop=12
num_loops=64
num_rx_per_device=4
num_devices=1

filepath = r"C:\ti\mmwave_studio_02_01_01_00\mmWaveStudio\PostProc\cascade_cpu_head_1\\*.bin"


def get_valid_num_frames(adc_idx_file_name):
    header_info_size = 6
    with open(adc_idx_file_name, 'rb') as idx_file:
        header_info = np.fromfile(idx_file, dtype=np.uint32, count=header_info_size)
        num_idx = header_info[3]  
    
    header_info_size = 3
    with open(adc_idx_file_name, 'rb') as idx_file:
        header_info = np.fromfile(idx_file, dtype=np.uint64, count=header_info_size)
        data_file_size = header_info[2]  
    
    return num_idx, data_file_size


data_files = glob.glob(filepath)
num_frames = 0 
data_file_size = 0
for file in data_files:
    if 'master_0000_idx' in file:
        num_frames, data_file_size = get_valid_num_frames(file)
        print(f"Number of Frames: {num_frames}")
        print(f"Data File Size: {data_file_size}")


angleFFTs_sum_frames = np.zeros((num_frames, 256, 256))
for frame_no in range(30):
    adc_data = read_adc_bin_tda2_separate_files(data_files,2, num_sample_per_chirp,num_chirp_per_loop,num_loops,num_rx_per_device,num_devices)
    rangeFFTOut = np.zeros(shape=adc_data.shape, dtype=np.complex128)
    DopplerFFTOut = np.zeros(shape=adc_data.shape, dtype=np.complex128)
    for i_tx in range(adc_data.shape[3]):
        rangeFFTOut[:,:,:,i_tx] = rangeProc(adc_data[:,:,:,i_tx])
        DopplerFFTOut[:,:,:,i_tx] = dopplerProc(rangeFFTOut[:,:,:,i_tx])
    DopplerFFTOut =  np.reshape(
        DopplerFFTOut, 
        (DopplerFFTOut.shape[0], DopplerFFTOut.shape[1], DopplerFFTOut.shape[2] * DopplerFFTOut.shape[3])
    )
    sig_integrate = 10 * np.log10(np.sum(np.abs(DopplerFFTOut)**2, axis=2) + 1)
    break


cfar_case(DopplerFFTOut)