import numpy as np
import pandas as pd
from obspy import read
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import os
from scipy.signal import butter, filtfilt


# Directory where seismic data files are stored
directory_path = r'C:\Users\Korisnik\Desktop\space_apps_2024_seismic_detection\space_apps_2024_seismic_detection\data\mars\test\data'
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

    def detect_and_filter_triggers(rms_energy, trigger_on, trigger_off):
        # Identify trigger points based on RMS energy
        on_off = []
        triggered = False

        for i in range(len(rms_energy)):
            if rms_energy[i] > trigger_on and not triggered:
                on_off.append((tr_times[i], 'on'))  # Log trigger ON
                triggered = True
            elif rms_energy[i] < trigger_off and triggered:
                on_off.append((tr_times[i], 'off'))  # Log trigger OFF
                triggered = False

        # Filter trigger points based on specified conditions
        filtered_on_off = []
        first_trigger_on = None
        last_trigger_off = None

        for j in range(len(on_off)):
            if on_off[j][1] == 'on':
                if first_trigger_on is None:
                    first_trigger_on = on_off[j]  # Store the first trigger ON
                if last_trigger_off is not None and (on_off[j][0] - last_trigger_off[0] > 100):
                    filtered_on_off.append(on_off[j])  # Append if more than 100 seconds from last OFF
            elif on_off[j][1] == 'off':
                last_trigger_off = on_off[j]  # Keep track of the last trigger OFF
                if (j + 1 < len(on_off)) and (on_off[j + 1][0] - last_trigger_off[0] > 100):
                    filtered_on_off.append(last_trigger_off)  # Append if more than 100 seconds from next ON

        # Always include the first trigger ON and the last trigger OFF
        if first_trigger_on:
            filtered_on_off.insert(0, first_trigger_on)  # Ensure the first ON is included
        if last_trigger_off:
            filtered_on_off.append(last_trigger_off)  # Ensure the last OFF is included

        # Add logic to filter triggers based on the time difference
        filtered_triggers = []
        i = 0
        while i < len(filtered_on_off):
            if i > 0:
                time_diff = filtered_on_off[i][0] - filtered_on_off[i - 1][0]
                # If the current trigger is OFF and the previous one is ON, check the time difference
                if (filtered_on_off[i][1] == 'off' and filtered_on_off[i - 1][1] == 'on' and time_diff < 40):
                    filtered_triggers.pop()  # Remove the last trigger ON
                    # Do not add the OFF trigger
                else:
                    filtered_triggers.append(filtered_on_off[i])
            else:
                filtered_triggers.append(filtered_on_off[i])
            i += 1

        return filtered_triggers

    # Apply the bandpass filter to the seismic data
    tr_data_filt = bandpass_filter(tr_data.astype(np.float64), minfreq, maxfreq, df)

    # Calculate RMS energy
    window_size = int(1 * df)  # 1-second window
    rms_energy = calculate_rms(tr_data_filt, window_size)

    # Check the range of RMS energy values
    print(f"RMS Energy values - Min: {np.min(rms_energy)}, Max: {np.max(rms_energy)}, Mean: {np.mean(rms_energy)}")

    # Set trigger thresholds with static values
    thr_on = 800  # Trigger on threshold
    thr_off = 400  # Trigger off threshold

    # Loop to handle static and dynamic trigger detection
    filtered_triggers = []
    for _ in range(2):  # Two passes: first for static, second for dynamic
        if len(filtered_triggers) == 0:  # Only detect if no triggers yet
            if _ == 0:  # Static thresholds
                current_thr_on = thr_on
                current_thr_off = thr_off
            else:  # Dynamic thresholds
                mean_rms = np.mean(rms_energy)
                std_rms = np.std(rms_energy)
                
                multiplier_on = 5
                multiplier_off = 2.5
                
                while len(filtered_triggers) == 0 and multiplier_on > 1 and multiplier_off > 1:
                    current_thr_on = mean_rms + (multiplier_on - 1) * std_rms
                    current_thr_off = mean_rms + (multiplier_off - 1) * std_rms
                    
                    print(f"Trying thresholds: ON={current_thr_on}, OFF={current_thr_off}")
                    
                    # Detect and filter triggers
                    filtered_triggers = detect_and_filter_triggers(rms_energy, current_thr_on, current_thr_off)
                    
                    # Decrease multipliers for next iteration
                    multiplier_on -= 1
                    multiplier_off -= 0.5  # Adjusting to reduce more gradually

                if len(filtered_triggers) == 0:
                    print("No triggers found with adjusted dynamic thresholds.")

            print(f"Using thresholds: ON={current_thr_on}, OFF={current_thr_off}")  # Debugging line

            # Detect and filter triggers
            filtered_triggers = detect_and_filter_triggers(rms_energy, current_thr_on, current_thr_off)

    # Create a figure with four subplots
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

    # Check if any triggers were detected in the second pass (dynamic detection)
    if len(filtered_triggers) == 0:
        mean_rms = np.mean(rms_energy)
        std_rms = np.std(rms_energy)
        thr_on_dynamic = mean_rms + (5 - 1) * std_rms
        thr_off_dynamic = mean_rms + (2.5 - 1) * std_rms
        axs[2].axhline(y=thr_on_dynamic, color='orange', linestyle='--', label='Trigger On Threshold (Dynamic)')
        axs[2].axhline(y=thr_off_dynamic, color='blue', linestyle='--', label='Trigger Off Threshold (Dynamic)')
    else:
        axs[2].axhline(y=current_thr_on, color='red', linestyle='--', label='Trigger On Threshold (Dynamic)')
        axs[2].axhline(y=current_thr_off, color='purple', linestyle='--', label='Trigger Off Threshold (Dynamic)')

    axs[2].set_title('RMS Energy with Trigger Thresholds')
    axs[2].set_xlabel('Time (s)')
    axs[2].set_ylabel('RMS Energy')
    axs[2].legend()

    # Plot 4: Raw Seismic Data with Current Trigger Markers
    axs[3].plot(tr_times, tr_data, label='Raw Seismic Data', alpha=0.5)
    for point in filtered_triggers:  # Only current triggers from the last detection
        if point[1] == 'on':
            axs[3].axvline(x=point[0], color='red', linestyle='--', label='Trigger ON' if 'Trigger ON' not in [line.get_label() for line in axs[3].get_lines()] else "")
        elif point[1] == 'off':
            axs[3].axvline(x=point[0], color='purple', linestyle='--', label='Trigger OFF' if 'Trigger OFF' not in [line.get_label() for line in axs[3].get_lines()] else "")
    axs[3].set_title('Raw Seismic Data with Current Trigger Markers')
    axs[3].set_xlabel('Time (s)')
    axs[3].set_ylabel('Amplitude')
    axs[3].legend()

    # Set x-ticks every 5000 seconds for all plots
    for ax in axs:
        ax.set_xticks(np.arange(min(tr_times), max(tr_times) + 1, 5000))
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.0f}'))  # Format ticks as integers

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
        a = 4
        b = 1
    
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
                magnitude = a+b*np.log10(abs(peak_amplitude*0.01))  # Just an example formula
                
                # Store the peak magnitude with corresponding times
                magnitude_list.append(magnitude)
                
                print(f"Peak magnitude: {magnitude}")
    
        return magnitude_list
    
    magnitude_list = find_peak_magnitude(tr_data, tr_times, filtered_triggers)
    print(magnitude_list)
    
    for mag in magnitude_list:
        print(f"Magnitude: {mag}")
            
    
        quake_data = {                      #default parameters which will be initialized
        "latitude": 4.502, 
        "longitude": 135.623, 
        "station": "The Insight Lander", 
        "magnitude": [4.5], 
        "date": "08/10/24", 
        "trigger_on": [], 
        "trigger_off": []
    }
        
    quake_data['magnitude'] = magnitude_list
    quake_data['date'] = stream[0].stats.starttime
    
    for t in filtered_triggers:
        if t[1] == 'on':  # If it's a trigger on point
            quake_data['trigger_on'].append(t[0])  # Append to trigger_on list
        elif t[1] == 'off':  # If it's a trigger off point
            quake_data['trigger_off'].append(t[0])  # Append to trigger_off list
    
    print(quake_data)


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
            save_directory = r'C:\Users\Korisnik\Desktop\space_apps_2024_seismic_detection\space_apps_2024_seismic_detection\data\mars\test\catalogs'
            catalog_filename = os.path.join(save_directory, 'test_data_catalogue_mars.csv')
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