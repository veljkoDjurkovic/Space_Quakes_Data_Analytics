import numpy as np
import pandas as pd
from obspy import read
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import os
from scipy.signal import butter, filtfilt


# Directory where seismic data files are stored
directory_path = r'C:\Users\Korisnik\Desktop\space_apps_2024_seismic_detection\space_apps_2024_seismic_detection\data\lunar\test\data\S16_GradeB'
name_list = {}

# Collect all the .mseed files from the directory
for filename in os.listdir(directory_path):
    if filename.endswith('.mseed'):
        file_path = os.path.join(directory_path, filename)
        name_list[filename] = file_path
        
for filename in name_list:
    
    try:
        mseed_file = name_list[filename]
        stream = read(mseed_file)
    except Exception as e:
        print(f"Error processing file {filename}: {e}")
        continue

    tr = stream.traces[0].copy() 
    print(tr.stats)
    tr_times = tr.times()  # Getting time values
    tr_data = tr.data
    starttime = tr.stats.starttime.datetime
    print('Successful iteration')

    # Define the sampling frequency of the trace
    df = tr.stats.sampling_rate

    # Bandpass filter parameters
    minfreq = 0.5  # Minimum frequency in Hz
    maxfreq = 1.0  # Maximum frequency in Hz

    # Bandpass filter function
    def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        b, a = butter(order, [low, high], btype='band')
        return b, a

    def bandpass_filter(data, lowcut, highcut, fs, order=5):
        b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        y = filtfilt(b, a, data)
        return y

    # Calculate RMS energy using a sliding window
    def calculate_rms(data, window_size):
        rms = np.array([np.sqrt(np.mean(np.square(data[i:i + window_size]))) 
                        for i in range(len(data) - window_size)])
        return rms

    # Detect triggers based on RMS energy and filter them
    def detect_and_filter_triggers(rms_energy, trigger_on, trigger_off):
        on_off = []
        triggered = False

        # Detect trigger points
        for i in range(len(rms_energy)):
            if rms_energy[i] > trigger_on and not triggered:
                on_off.append((tr_times[i], 'on'))  # Log trigger ON
                triggered = True
            elif rms_energy[i] < trigger_off and triggered:
                on_off.append((tr_times[i], 'off'))  # Log trigger OFF
                triggered = False

        # Filter triggers based on time gaps
        filtered_on_off = []
        first_trigger_on = None
        last_trigger_off = None

        for j in range(len(on_off)):
            if on_off[j][1] == 'on':
                if first_trigger_on is None:
                    first_trigger_on = on_off[j]  # Store the first trigger ON
                if last_trigger_off is not None and (on_off[j][0] - last_trigger_off[0] > 700):
                    filtered_on_off.append(on_off[j])  # Append if more than 700 seconds from last OFF
            elif on_off[j][1] == 'off':
                last_trigger_off = on_off[j]  # Keep track of the last trigger OFF
                if (j + 1 < len(on_off)) and (on_off[j + 1][0] - last_trigger_off[0] > 700):
                    filtered_on_off.append(last_trigger_off)  # Append if more than 700 seconds from next ON

        # Always include the first ON and last OFF triggers
        if first_trigger_on:
            filtered_on_off.insert(0, first_trigger_on)
        if last_trigger_off:
            filtered_on_off.append(last_trigger_off)

        # Filter triggers based on time difference
        filtered_triggers = []
        i = 0
        while i < len(filtered_on_off):
            if i > 0:
                time_diff = filtered_on_off[i][0] - filtered_on_off[i - 1][0]
                if (filtered_on_off[i][1] == 'off' and filtered_on_off[i - 1][1] == 'on' and time_diff < 300):
                    filtered_triggers.pop()  # Remove last trigger ON
                else:
                    filtered_triggers.append(filtered_on_off[i])
            else:
                filtered_triggers.append(filtered_on_off[i])
            i += 1

        return filtered_triggers

    # Apply bandpass filter to the seismic data
    tr_data_filt = bandpass_filter(tr_data.astype(np.float64), minfreq, maxfreq, df)

    # Calculate RMS energy
    window_size = int(1 * df)  # 1-second window
    rms_energy = calculate_rms(tr_data_filt, window_size)

    # Print range of RMS energy values for debugging
    print(f"RMS Energy values - Min: {np.min(rms_energy)}, Max: {np.max(rms_energy)}, Mean: {np.mean(rms_energy)}")

    # Set static trigger thresholds
    thr_on = 1.5e-9
    thr_off = 1.0e-9

    # First detect static and dynamic triggers
    filtered_triggers = []
    for _ in range(2):  # Two passes: first static, then dynamic
        if len(filtered_triggers) == 0:  # Detect if no triggers yet
            if _ == 0:  # Static detection
                current_thr_on = thr_on
                current_thr_off = thr_off
            else:  # Dynamic detection based on RMS statistics
                mean_rms = np.mean(rms_energy)
                std_rms = np.std(rms_energy)
                current_thr_on = mean_rms + 5 * std_rms
                current_thr_off = mean_rms + 2.5 * std_rms

            print(f"Using thresholds: ON={current_thr_on}, OFF={current_thr_off}")

            # Detect and filter triggers
            filtered_triggers = detect_and_filter_triggers(rms_energy, current_thr_on, current_thr_off)

    # Plotting the results in four subplots
    fig, axs = plt.subplots(4, 1, figsize=(12, 24))

    # Plot 1: Raw Seismic Data
    axs[0].plot(tr_times, tr_data, label='Raw Seismic Data', alpha=0.5)
    axs[0].set_title('Raw Seismic Data')
    axs[0].set_xlabel('Time (s)')
    axs[0].set_ylabel('Amplitude')
    axs[0].legend()

    # Plot 2: Filtered Seismic Data
    axs[1].plot(tr_times, tr_data_filt, label='Filtered Seismic Data', color='orange', alpha=0.5)
    axs[1].set_title('Filtered Seismic Data')
    axs[1].set_xlabel('Time (s)')
    axs[1].set_ylabel('Amplitude')
    axs[1].legend()

    # Plot 3: RMS Energy with Trigger Thresholds
    axs[2].plot(tr_times[:len(rms_energy)], rms_energy, label='RMS Energy', color='green')

    # Dynamic threshold if no triggers detected
    if len(filtered_triggers) == 0:
        mean_rms = np.mean(rms_energy)
        std_rms = np.std(rms_energy)
        thr_on_dynamic = mean_rms + 5 * std_rms
        thr_off_dynamic = mean_rms + 2.5 * std_rms
        axs[2].axhline(y=thr_on_dynamic, color='orange', linestyle='--', label='Trigger On Threshold (Dynamic)')
        axs[2].axhline(y=thr_off_dynamic, color='blue', linestyle='--', label='Trigger Off Threshold (Dynamic)')
    else:
        axs[2].axhline(y=current_thr_on, color='red', linestyle='--', label='Trigger On Threshold')
        axs[2].axhline(y=current_thr_off, color='purple', linestyle='--', label='Trigger Off Threshold')

    axs[2].set_title('RMS Energy with Trigger Thresholds')
    axs[2].set_xlabel('Time (s)')
    axs[2].set_ylabel('RMS Energy')
    axs[2].legend()

    # Plot 4: Raw Seismic Data with Trigger Markers
    axs[3].plot(tr_times, tr_data, label='Raw Seismic Data', alpha=0.5)
    for point in filtered_triggers:
        if point[1] == 'on':
            axs[3].axvline(x=point[0], color='red', linestyle='--', label='Trigger ON' if 'Trigger ON' not in [line.get_label() for line in axs[3].get_lines()] else "")
        elif point[1] == 'off':
            axs[3].axvline(x=point[0], color='blue', linestyle='--', label='Trigger OFF' if 'Trigger OFF' not in [line.get_label() for line in axs[3].get_lines()] else "")
    
    axs[3].set_title('Raw Seismic Data with Trigger Points')
    axs[3].set_xlabel('Time (s)')
    axs[3].set_ylabel('Amplitude')
    axs[3].legend()

   # Set x-ticks every 5000 seconds for all plots
    for ax in axs:
        ax.set_xticks(np.arange(min(tr_times), max(tr_times) + 1, 5000))

    plt.tight_layout()
    plt.show()

    # Print plotted triggers
    print("\nPlotted Triggers:")
    if not filtered_triggers:
        print("No filtered triggers found.")
    else:
        for point in filtered_triggers:
            print(f"{point[1].capitalize()} at: {point[0]}")
    
    def find_peak_magnitude(trace_data, trace_times, filtered_triggers):
    # Iterate through triggers to find 'on' and 'off' points
        magnitude_list = []
        a = 6.13
        b = 2/3
    
        for i in range(1, len(filtered_triggers)):
            if filtered_triggers[i-1][1] == 'on' and filtered_triggers[i][1] == 'off':
                # Extract indices corresponding to 'on' and 'off' times
                on_time = filtered_triggers[i-1][0]
                off_time = filtered_triggers[i][0]
                
                # Find the indices in the trace_times that match these trigger points
                on_idx = np.searchsorted(trace_times, on_time)
                off_idx = np.searchsorted(trace_times, off_time)
                
                # Slice the data between the trigger on and off points
                sliced_data = trace_data[on_idx:off_idx]
                sliced_times = trace_times[on_idx:off_idx]
                
                # Find the maximum amplitude (peak value) in the sliced data
                peak_amplitude = np.max(np.abs(sliced_data))
                
                # Estimate the magnitude (simple example, you would need a more specific relation for real data)
                magnitude = a+b*np.log10(abs(peak_amplitude))  # Just an example formula
                
                # Store the peak magnitude with corresponding times
                magnitude_list.append(magnitude)
                
                print(f"Peak magnitude: {magnitude}")
    
        return magnitude_list
    
    magnitude_list = find_peak_magnitude(tr_data, tr_times, filtered_triggers)
    print(magnitude_list)
    
    for mag in magnitude_list:
        print(f"Magnitude: {mag}")
        
    #Making dictionary for MySQL Data base
           
  
    quake_data = {                      #default parameters which will be initialized
        "latitude": 0.674, 
        "longitude": 23.473, 
        "station": "Apollo 11 - Mare Tranquillitatis", 
        "magnitude": [4.5], 
        "date": "08/10/24", 
        "trigger_on": [], 
        "trigger_off": []
    }
    
    if stream[0].stats.station == 'S11':
        quake_data['latitude'] = 0.674
        quake_data['longitude'] = 23.473
        quake_data['station'] = 'S11'
        
    elif stream[0].stats.station == 'S12':
        quake_data['latitude'] = -3.012
        quake_data['longitude'] = -34.434
        quake_data['station'] = 'S12'
    
    elif stream[0].stats.station == 'S14':  
        quake_data['latitude'] = -4.1
        quake_data['longitude'] = -17.5
        quake_data['station'] = 'S14'
        
    elif stream[0].stats.station == 'S15': 
        quake_data['latitude'] = -26.1
        quake_data['longitude'] = 3.6
        quake_data['station'] = 'S15'
    
    elif stream[0].stats.station == 'S16': 
        quake_data['latitude'] = -8.97
        quake_data['longitude'] = 15.5
        quake_data['station'] = 'S16'
    
    elif stream[0].stats.station == 'S17':  
        quake_data['latitude'] = 20.1
        quake_data['longitude'] = 30.8
        quake_data['station'] = 'S17'
    else:
        quake_data['latitude'] = 0
        quake_data['longitude'] = 0
        quake_data['station'] = 'Undefined'
        
    quake_data['magnitude'] = magnitude_list
    quake_data['date'] = stream[0].stats.starttime
    
    for t in filtered_triggers:
        if t[1] == 'on':  # If it's a trigger on point
            quake_data['trigger_on'].append(t[0])  # Append to trigger_on list
        elif t[1] == 'off':  # If it's a trigger off point
            quake_data['trigger_off'].append(t[0])  # Append to trigger_off list
    
    print(quake_data)
    
    
    # Function to write detected events into a catalog (DataFrame)
    def write_catalog(trace, filtered_triggers, trace_times, filename):
        # Get the start time of the trace
        starttime = trace.stats.starttime.datetime
        
        # Iterate through triggers to compile detection times
        detection_times = []
        fnames = []
        relative_times = []

        for i in range(1, len(filtered_triggers)):
            if filtered_triggers[i-1][1] == 'on' and filtered_triggers[i][1] == 'off':
                on_time = starttime + timedelta(seconds=trace_times[np.searchsorted(trace_times, filtered_triggers[i-1][0])])
                on_time_str = datetime.strftime(on_time, '%Y-%m-%dT%H:%M:%S.%f')
                
                # Append data for catalog
                detection_times.append(on_time_str)
                fnames.append(filename)
                relative_times.append(trace_times[np.searchsorted(trace_times, filtered_triggers[i-1][0])])

        # Compile into a DataFrame
        catalog_df = pd.DataFrame({
            'filename': fnames,
            'time_abs(%Y-%m-%dT%H:%M:%S.%f)': detection_times,
            'time_rel(sec)': relative_times
        })
    
        # Specify the directory where CSV file will be saved
        save_directory = r'C:\Users\Korisnik\Desktop\space_apps_2024_seismic_detection\space_apps_2024_seismic_detection\data\lunar\test\catalogs'
        catalog_filename = os.path.join(save_directory, 's16_gradeB_catalog_moon.csv')
        # Create the full path for the CSV file
        # Check if the file already exists to handle header
        if not os.path.exists(catalog_filename):
            # If file doesn't exist, write header and data
            catalog_df.to_csv(catalog_filename, index=False, mode='w', header=True)
        else:
            # If file exists, append data without writing header
            catalog_df.to_csv(catalog_filename, index=False, mode='a', header=False)
            
        print(f"Catalog saved to {catalog_filename}")
        return catalog_df
        
    # write_catalog(tr, filtered_triggers, tr_times, filename)

    