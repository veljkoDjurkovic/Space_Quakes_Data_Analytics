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

def insert_earthquake_data(data):
    try:
        connection = mysql.connector.connect(**db_config)

        if connection.is_connected():
            cursor = connection.cursor()
            insert_query = """
            INSERT INTO nasa (latitude, longitude, station, magnitude, date, trigger_on, trigger_off)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
            """
            for entry in data:
                cursor.execute(insert_query, (
                    entry["latitude"], 
                    entry["longitude"], 
                    entry["station"], 
                    entry["magnitude"], 
                    entry["date"], 
                    entry["trigger_on"], 
                    entry["trigger_off"]
                ))
            
            connection.commit()  
            print("Earthquake data inserted successfully")

    except Error as e:
        print(f"Error while connecting to MySQL: {e}")
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()

@app.route('/earthquake-data', methods=['GET'])
def get_earthquake_data():
    try:
        connection = mysql.connector.connect(**db_config)
        if connection.is_connected():
            cursor = connection.cursor(dictionary=True)
            cursor.execute("SELECT * FROM moon")
            earthquake_data = cursor.fetchall()
            return jsonify(earthquake_data)
    except Error as e:
        return jsonify({"error": str(e)})
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()

@app.route('/mars-earthquake-data', methods=['GET'])
def get_mars_earthquake_data():
    try:
        mars_db_config = db_config = {
        'host': 'localhost',  
        'database': 'nasa',  
        'user': 'root',  
        'password': '' 
        }

        connection = mysql.connector.connect(**mars_db_config)
        if connection.is_connected():
            cursor = connection.cursor(dictionary=True)
            cursor.execute("SELECT * FROM mars")
            mars_earthquake_data = cursor.fetchall()
            return jsonify(mars_earthquake_data)
    except Error as e:
        return jsonify({"error": str(e)})
    finally:
        if connection.is_connected():
            cursor.close()
            connection.close()

if __name__ == '__main__':
    app.run(debug=True)
