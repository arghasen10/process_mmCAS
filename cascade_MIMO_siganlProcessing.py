import glob
import numpy as np
from utils.rangeProcCascade import rangeProc
from utils.DopplerProcCascade import dopplerProc
import matplotlib.pyplot as plt
import seaborn as sns

num_sample_per_chirp=256
num_chirp_per_loop=12
num_loops=64
num_rx_per_device=4
num_devices=1
chunk_size = num_sample_per_chirp*num_chirp_per_loop*num_loops*num_rx_per_device*num_devices*4
RxForMimoProcess = [12,13,14,15,0,1,2,3,8,9,10,11,4,5,6,7]

def bin2np_frame(bin_frame):  #
    neg = np.bitwise_and(bin_frame, 1 << 15) != 0
    bin_frame[neg] = bin_frame[neg]-(2**16)
    bin_frame = bin_frame[0::2] + 1j * bin_frame[1::2]
    bin_frame = np.reshape(bin_frame, (num_rx_per_device, num_sample_per_chirp, num_chirp_per_loop, num_loops), order='F')
    bin_frame = np.transpose(bin_frame,axes=(1,3,0,2))
    return bin_frame

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
file_paths = {}
dtype=np.uint16
for file in data_files:
    if 'data' in file:
        if 'master' in file:
            file_paths['master'] = file
        elif 'slave1' in file:
            file_paths['slave1'] = file
        elif 'slave2' in file:
            file_paths['slave2'] = file
        elif 'slave3' in file:
            file_paths['slave3'] = file
with open(file_paths['master'], 'rb') as f_master, \
     open(file_paths['slave1'], 'rb') as f_slave1, \
     open(file_paths['slave2'], 'rb') as f_slave2, \
     open(file_paths['slave3'], 'rb') as f_slave3:
    for frame_no in range(30):
        masterArr = bin2np_frame(np.frombuffer(f_master.read(chunk_size), dtype=np.int16).copy())
        slave1_arr = bin2np_frame(np.frombuffer(f_slave1.read(chunk_size), dtype=np.int16).copy())
        slave2_arr = bin2np_frame(np.frombuffer(f_slave2.read(chunk_size), dtype=np.int16).copy())
        slave3_arr = bin2np_frame(np.frombuffer(f_slave3.read(chunk_size), dtype=np.int16).copy())
        radar_cube=np.zeros(shape=(256,64,16,12), dtype=np.complex_)
        radar_cube[:,:,0:4,:] = masterArr
        radar_cube[:,:,4:8,:] = slave1_arr
        radar_cube[:,:,8:12,:] = slave2_arr
        radar_cube[:,:,12:,:] = slave3_arr
        radar_cube = radar_cube[:,:,RxForMimoProcess,:]
        print(radar_cube.shape)
        print(radar_cube[0,0,1,0:10])
        rangeFFTOut = np.zeros(shape=radar_cube.shape, dtype=np.complex_)
        DopplerFFTOut = np.zeros(shape=radar_cube.shape, dtype=np.complex_)
        for i_tx in range(radar_cube.shape[3]):
            rangeFFTOut[:,:,:,i_tx] = rangeProc(radar_cube[:,:,:,i_tx])
            DopplerFFTOut[:,:,:,i_tx] = dopplerProc(rangeFFTOut[:,:,:,i_tx])
        DopplerFFTOut =  np.reshape(
            DopplerFFTOut, 
            (DopplerFFTOut.shape[0], DopplerFFTOut.shape[1], DopplerFFTOut.shape[2] * DopplerFFTOut.shape[3])
        )
        sig_integrate = 10 * np.log10(np.sum(np.abs(DopplerFFTOut)**2, axis=2) + 1)
        if frame_no == 2:
            break


sns.heatmap(np.abs(np.sum(DopplerFFTOut,axis=(2))))
plt.show()
print(np.abs(np.sum(DopplerFFTOut,axis=(2)))[:10,0])
 
