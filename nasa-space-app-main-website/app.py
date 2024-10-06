from flask import Flask, jsonify, request
from flask_cors import CORS
import mysql.connector
from mysql.connector import Error

import numpy as np
import pandas as pd
from obspy import read
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import os
from scipy.signal import butter, filtfilt

app = Flask(__name__)
CORS(app)

db_config = {
    'host': 'localhost',  
    'database': 'nasa',  
    'user': 'root',  
    'password': '' 
}

def insert_earthquake_data(data):
    try:
        connection = mysql.connector.connect(**db_config)

        if connection.is_connected():
            cursor = connection.cursor()
            insert_query = """
            INSERT INTO moon (latitude, longitude, station, magnitude, date, trigger_on, trigger_off)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
            """
            for entry in data:
                cursor.execute(insert_query, (
                    entry['latitude'], 
                    entry['longitude'], 
                    entry['station'], 
                    entry['magnitude'], 
                    entry['date'], 
                    entry['trigger_on'], 
                    entry['trigger_off']
                ))
            
            connection.commit()  
            print("Earthquake data inserted successfully")

    except Error as e:
        print(f"Error while connecting to MySQL: {e}")
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()



def gather_and_insert_earthquake_data():

    directory_path = r'C:\Users\Korisnik\Desktop\space_apps_2024_seismic_detection\space_apps_2024_seismic_detection\data\lunar\test\data\S16_GradeA'

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

        trace = stream.traces[0].copy()  
        print(trace.stats)
        trace_times = trace.times() 
        trace_data = trace.data
        starttime = trace.stats.starttime.datetime
        print('Successful iteration')

        df = trace.stats.sampling_rate

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
                    on_off.append((trace_times[i], 'on'))  
                    triggered = True
                elif rms_energy[i] < trigger_off and triggered:
                    on_off.append((trace_times[i], 'off'))  
                    triggered = False

            filtered_on_off = []
            first_trigger_on = None
            last_trigger_off = None

            for j in range(len(on_off)):
                if on_off[j][1] == 'on':
                    if first_trigger_on is None:
                        first_trigger_on = on_off[j]  
                    if last_trigger_off is not None and (on_off[j][0] - last_trigger_off[0] > 700):
                        filtered_on_off.append(on_off[j])  
                elif on_off[j][1] == 'off':
                    last_trigger_off = on_off[j]  
                    if (j + 1 < len(on_off)) and (on_off[j + 1][0] - last_trigger_off[0] > 700):
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
                    if (filtered_on_off[i][1] == 'off' and filtered_on_off[i - 1][1] == 'on' and time_diff < 300):
                        filtered_triggers.pop()  
                    else:
                        filtered_triggers.append(filtered_on_off[i])
                else:
                    filtered_triggers.append(filtered_on_off[i])
                i += 1

            return filtered_triggers

        tr_data_filt = bandpass_filter(trace_data.astype(np.float64), minfreq, maxfreq, df)

        window_size = int(1 * df) 
        rms_energy = calculate_rms(tr_data_filt, window_size)

        print(f"RMS Energy values - Min: {np.min(rms_energy)}, Max: {np.max(rms_energy)}, Mean: {np.mean(rms_energy)}")

        thr_on = 1.5e-9
        thr_off = 1.0e-9

        filtered_triggers = []
        for _ in range(2): 
            if len(filtered_triggers) == 0:  
                if _ == 0:  
                    current_thr_on = thr_on
                    current_thr_off = thr_off
                else:  
                    mean_rms = np.mean(rms_energy)
                    std_rms = np.std(rms_energy)
                    current_thr_on = mean_rms + 5 * std_rms
                    current_thr_off = mean_rms + 2.5 * std_rms

                print(f"Using thresholds: ON={current_thr_on}, OFF={current_thr_off}")

                filtered_triggers = detect_and_filter_triggers(rms_energy, current_thr_on, current_thr_off)


        def find_peak_magnitude(trace_data, trace_times, filtered_triggers):
            magnitude_list = []
            a = 6.13
            b = 0.66
        
            for i in range(1, len(filtered_triggers)):
                if filtered_triggers[i-1][1] == 'on' and filtered_triggers[i][1] == 'off':
                    on_time = filtered_triggers[i-1][0]
                    off_time = filtered_triggers[i][0]
                    
                    on_idx = np.searchsorted(trace_times, on_time)
                    off_idx = np.searchsorted(trace_times, off_time)
                    
                    sliced_data = trace_data[on_idx:off_idx]
                    peak_amplitude = np.max(np.abs(sliced_data))
                    
                    magnitude = a + b * np.log10(abs(peak_amplitude))
                    magnitude_list.append(magnitude)
                    print(f"Peak magnitude: {magnitude}")

            return magnitude_list
        
        magnitude_list = find_peak_magnitude(trace_data, trace_times, filtered_triggers)
        print(magnitude_list)

        quake_data = {                      
            "latitude": 0.674, 
            "longitude": 23.473, 
            "station": "Apollo 11 - Mare Tranquillitatis", 
            "magnitude": None,  
            "date": stream[0].stats.starttime.strftime("%Y-%m-%d %H:%M:%S"), 
            "trigger_on": None, 
            "trigger_off": None
        }

        if stream[0].stats.station == 'S11':
            quake_data["latitude"] = 0.674
            quake_data["longitude"] = 23.473
            quake_data["station"] = "Apollo 11 - Mare Tranquillitatis"
        elif stream[0].stats.station == 'S12':
            quake_data["latitude"] = -3.012
            quake_data["longitude"] = -34.434
            quake_data["station"] = "Apollo 12 - Oceanus Procellarum"
        elif stream[0].stats.station == 'S14':  
            quake_data["latitude"] = -4.1
            quake_data["longitude"] = -17.5
            quake_data["station"] = "Apollo 14 - Fra Mauro Highlands"
        elif stream[0].stats.station == 'S15': 
            quake_data["latitude"] = -26.1
            quake_data["longitude"] = 3.6
            quake_data["station"] = "Apollo 15 - Hadley-Apennine"
        elif stream[0].stats.station == 'S16': 
            quake_data["latitude"] = -8.97
            quake_data["longitude"] = 15.5
            quake_data["station"] = "Apollo 16 - Descartes Highlands"
        elif stream[0].stats.station == 'S17':  
            quake_data["latitude"] = -8.2
            quake_data["longitude"] = 22.3
            quake_data["station"] = "Apollo 17 - Taurus-Littrow"

        if magnitude_list:  
            quake_data["magnitude"] = float(magnitude_list[0])  
        else:
            quake_data["magnitude"] = None 

        if filtered_triggers:
            quake_data["trigger_on"] = float(filtered_triggers[0][0])  
            quake_data["trigger_off"] = float(filtered_triggers[-1][0]) 
        else:
            quake_data["trigger_on"] = None 
            quake_data["trigger_off"] = None

        print(quake_data)

        insert_earthquake_data([quake_data])  

    
    # insert_earthquake_data(earthquake_data) 

@app.route('/earthquake-data', methods=['GET'])
def get_earthquake_data():
    try:
        connection = mysql.connector.connect(**db_config)
        if connection.is_connected():
            cursor = connection.cursor(dictionary=True)
            cursor.execute("SELECT * FROM nasa")
            earthquake_data = cursor.fetchall()
            return jsonify(earthquake_data)
    except Error as e:
        return jsonify({"error": str(e)})
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()

if __name__ == '__main__':
    gather_and_insert_earthquake_data()  
    app.run(debug=True)
