import numpy as np
import glob

TxToEnable = [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
RxForMimoProcess = [12,13,14,15,0,1,2,3,8,9,10,11,4,5,6,7]
def read_bin_file(file_full_path, frame_idx, 
                  num_sample_per_chirp, num_chirp_per_loop, 
                  num_loops, num_rx_per_device, num_devices):
    expected_num_samples_per_frame = (num_sample_per_chirp * 
                                      num_chirp_per_loop * 
                                      num_loops * 
                                      num_rx_per_device * 2)
    
    with open(file_full_path, 'rb') as fp:
        # Seek to the frame position
        fp.seek((frame_idx - 1) * expected_num_samples_per_frame * 2, 0)
        
        # Read ADC data
        adc_data = np.fromfile(fp, dtype=np.uint16, count=expected_num_samples_per_frame)
        
    # Convert unsigned to signed (handle two's complement)
    neg = (adc_data & (1 << 15)) != 0  # Check if 16th bit is set
    adc_data[neg] -= (1 << 16)
    
    # Create complex numbers from interleaved real and imaginary parts
    adc_data = adc_data[0::2] + 1j * adc_data[1::2]
    
    # Reshape and reorder dimensions
    adc_data_complex = adc_data.reshape(num_rx_per_device, 
                                        num_sample_per_chirp, 
                                        num_chirp_per_loop, 
                                        num_loops)
    
    # Permute dimensions to match MATLAB order
    adc_data_complex = np.transpose(adc_data_complex, (1, 3, 0, 2))
    
    return adc_data_complex


def read_adc_bin_tda2_separate_files(data_files, frame_idx, num_sample_per_chirp,
                                     num_chirp_per_loop, num_loops, num_rx_per_device, num_devices):
    file_paths = {}
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

    
    radar_data = np.zeros((num_sample_per_chirp, num_loops, num_rx_per_device * 4, num_chirp_per_loop), dtype=np.complex64)
    
    radar_data[:, :, 0:4, :] = read_bin_file(file_paths['master'], frame_idx, num_sample_per_chirp, num_chirp_per_loop, num_loops, num_rx_per_device, num_devices)
    radar_data[:, :, 4:8, :] = read_bin_file(file_paths['slave1'], frame_idx, num_sample_per_chirp, num_chirp_per_loop, num_loops, num_rx_per_device, num_devices)
    radar_data[:, :, 8:12, :] = read_bin_file(file_paths['slave2'], frame_idx, num_sample_per_chirp, num_chirp_per_loop, num_loops, num_rx_per_device, num_devices)
    radar_data[:, :, 12:16, :] = read_bin_file(file_paths['slave3'], frame_idx, num_sample_per_chirp, num_chirp_per_loop, num_loops, num_rx_per_device, num_devices)
    radar_data = radar_data[:, :, :, TxToEnable]
    radar_data = radar_data[:,:,RxForMimoProcess,:]
    return radar_data

if __name__=="__main__":
    # file_full_path = r"C:\\ti\\mmwave_studio_02_01_01_00\\mmWaveStudio\\PostProc\\cascade_cpu_head_1\\master_0000_data.bin"
    # adc_data_complex = read_bin_file(file_full_path, 
    #                                 frame_idx=1, 
    #                                 num_sample_per_chirp=256, 
    #                                 num_chirp_per_loop=12, 
    #                                 num_loops=64, 
    #                                 num_rx_per_device=4, 
    #                                 num_devices=1)
    filepath = r"C:\ti\mmwave_studio_02_01_01_00\mmWaveStudio\PostProc\cascade_cpu_head_1\\*.bin"
    data_files = glob.glob(filepath)
    radar_data = read_adc_bin_tda2_separate_files(data_files,
                                     frame_idx=1, 
                                     num_sample_per_chirp=256,
                                     num_chirp_per_loop=12,
                                     num_loops=64,
                                     num_rx_per_device=4,
                                     num_devices=1)
    print(radar_data.shape)