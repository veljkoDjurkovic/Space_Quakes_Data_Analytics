import numpy as np
import pandas as pd
from obspy import read
from datetime import datetime, timedelta
import os
from scipy.signal import butter, filtfilt
from flask import Flask, jsonify, request
from flask_cors import CORS
import mysql.connector
from mysql.connector import Error

app = Flask(__name__)
CORS(app)

db_config = {
    'host': 'localhost', 
    'database': 'nasa', 
    'user': 'root',  
    'password': ''  
}

def insert_earthquake_data(entry):
    try:
        connection = mysql.connector.connect(**db_config)

        if connection.is_connected():
            cursor = connection.cursor()
            insert_query = """
            INSERT INTO mars (latitude, longitude, station, magnitude, date, trigger_on, trigger_off)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
            """
            cursor.execute(insert_query, (entry['latitude'], entry['longitude'], entry['station'], entry['magnitude'], entry['date'], entry['trigger_on'], entry['trigger_off']))
            
            connection.commit()  
            print("Earthquake data inserted successfully")

    except Error as e:
        print(f"Error while connecting to MySQL: {e}")
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()

def find_peak_magnitude(trace_data, trace_times, filtered_triggers):
    magnitude_list = []
    for i in range(1, len(filtered_triggers)):
        if filtered_triggers[i-1][1] == 'on' and filtered_triggers[i][1] == 'off':
            on_time = filtered_triggers[i-1][0]
            off_time = filtered_triggers[i][0]
            
            on_idx = np.searchsorted(trace_times, on_time)
            off_idx = np.searchsorted(trace_times, off_time)
            
            sliced_data = trace_data[on_idx:off_idx]
            peak_amplitude = np.max(np.abs(sliced_data))
            magnitude = 4.0 + 1.0 * np.log10(abs(peak_amplitude * 0.01))  
            magnitude_list.append(magnitude)
    return magnitude_list

def process_and_insert_data(tr_data, tr_times, filtered_triggers):
    magnitude_list = find_peak_magnitude(tr_data, tr_times, filtered_triggers)
    for i in range(len(magnitude_list)):
        quake_data = {
            "latitude": 4.502, 
            "longitude": 135.623, 
            "station": "The Insight Lander", 
            "date": datetime.strftime(datetime.now(), '%Y/%m/%d'),  
            "magnitude": float(magnitude_list[i]),  
            "trigger_on": float(filtered_triggers[i * 2][0]) if i * 2 < len(filtered_triggers) else None,  
            "trigger_off": float(filtered_triggers[i * 2 + 1][0]) if i * 2 + 1 < len(filtered_triggers) else None,  
        }
        insert_earthquake_data(quake_data)  

directory_path = r'C:\Users\Korisnik\Desktop\space_apps_2024_seismic_detection\space_apps_2024_seismic_detection\data\mars\test\data'
name_list = {}

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
    tr_times = tr.times()  
    tr_data = tr.data
    starttime = tr.stats.starttime.datetime

   
    df = tr.stats.sampling_rate

    minfreq = 0.5  
    maxfreq = 1.0  

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

    def calculate_rms(data, window_size):
        rms = np.array([np.sqrt(np.mean(np.square(data[i:i + window_size]))) 
                        for i in range(len(data) - window_size)])
        return rms

    def detect_and_filter_triggers(rms_energy, trigger_on, trigger_off):
        on_off = []
        triggered = False

        for i in range(len(rms_energy)):
            if rms_energy[i] > trigger_on and not triggered:
                on_off.append((tr_times[i], 'on')) 
                triggered = True
            elif rms_energy[i] < trigger_off and triggered:
                on_off.append((tr_times[i], 'off'))  
                triggered = False

        filtered_on_off = []
        first_trigger_on = None
        last_trigger_off = None

        for j in range(len(on_off)):
            if on_off[j][1] == 'on':
                if first_trigger_on is None:
                    first_trigger_on = on_off[j]  
                if last_trigger_off is not None and (on_off[j][0] - last_trigger_off[0] > 100):
                    filtered_on_off.append(on_off[j]) 
            elif on_off[j][1] == 'off':
                last_trigger_off = on_off[j]  
                if (j + 1 < len(on_off)) and (on_off[j + 1][0] - last_trigger_off[0] > 100):
                    filtered_on_off.append(last_trigger_off)  

        if first_trigger_on:
            filtered_on_off.insert(0, first_trigger_on)  
        if last_trigger_off:
            filtered_on_off.append(last_trigger_off)  

        filtered_triggers = []
        i = 0
        while i < len(filtered_on_off):
            if i > 0:
                time_diff = filtered_on_off[i][0] - filtered_on_off[i - 1][0]
                if (filtered_on_off[i][1] == 'off' and filtered_on_off[i - 1][1] == 'on' and time_diff < 40):
                    filtered_triggers.pop() 
                else:
                    filtered_triggers.append(filtered_on_off[i])
            else:
                filtered_triggers.append(filtered_on_off[i])
            i += 1

        return filtered_triggers

    tr_data_filt = bandpass_filter(tr_data.astype(np.float64), minfreq, maxfreq, df)

    window_size = int(1 * df)  
    rms_energy = calculate_rms(tr_data_filt, window_size)

    thr_on = 800  
    thr_off = 400  

    filtered_triggers = []
    for _ in range(2): 
        if len(filtered_triggers) == 0: 
            if _ == 0:  
                current_thr_on = thr_on
                current_thr_off = thr_off
            else:  
                mean_rms = np.mean(rms_energy)
                std_rms = np.std(rms_energy)
                
                multiplier_on = 5
                multiplier_off = 2.5
                
                while len(filtered_triggers) == 0 and multiplier_on > 1 and multiplier_off > 1:
                    current_thr_on = mean_rms + (multiplier_on - 1) * std_rms
                    current_thr_off = mean_rms + (multiplier_off - 1) * std_rms
                    
                    filtered_triggers = detect_and_filter_triggers(rms_energy, current_thr_on, current_thr_off)
                    
                    multiplier_on -= 1
                    multiplier_off -= 0.5 

                if len(filtered_triggers) == 0:
                    print("No triggers found with adjusted dynamic thresholds.")

            filtered_triggers = detect_and_filter_triggers(rms_energy, current_thr_on, current_thr_off)

    process_and_insert_data(tr_data, tr_times, filtered_triggers)

if __name__ == '__main__':
    app.run(debug=True)
