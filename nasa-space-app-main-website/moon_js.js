

function initMap() {
    const map = new google.maps.Map(document.getElementById("map"), {
        center: { lat: 0, lng: 0 },
        zoom: 2,
        streetViewControl: false,
        mapTypeControlOptions: {
            mapTypeIds: ["moon"],
        },
    });

    const moonMapType = new google.maps.ImageMapType({
        getTileUrl: function (coord, zoom) {
            const normalizedCoord = getNormalizedCoord(coord, zoom);
            if (!normalizedCoord) return "";
            const bound = Math.pow(2, zoom);
            return "https://mw1.google.com/mw-planetary/lunar/lunarmaps_v1/clem_bw" + "/" + zoom + "/" + normalizedCoord.x + "/" + (bound - normalizedCoord.y - 1) + ".jpg";
        },
        tileSize: new google.maps.Size(256, 256),
        maxZoom: 9,
        minZoom: 0,
        radius: 1738000,
        name: "Moon",
    });

    map.mapTypes.set("moon", moonMapType);
    map.setMapTypeId("moon");

    const stations = [
        { "lat": 0.674, "lng": 23.473, "location": "Apollo 11 - Mare Tranquillitatis", "ind": "S11" },
        { "lat": -3.012, "lng": -34.434, "location": "Apollo 12 - Oceanus Procellarum" , "ind": "S12"},
        { "lat": -4.1, "lng": -17.5, "location": "Apollo 14 - Fra Mauro Highlands" , "ind": "S14" },
        { "lat": 26.1, "lng": 3.6, "location": "Apollo 15 - Hadley-Apennine" , "ind": "S15"},
        { "lat": -8.97, "lng": 15.5, "location": "Apollo 16 - Descartes Highlands", "ind": "S16" },
        { "lat": 20.1, "lng": 30.8, "location": "Apollo 17 - Taurus-Littrow Valley" , "ind": "S17"}
    ];
    

    fetch('http://127.0.0.1:5000/earthquake-data')
    .then(response => response.json())
    .then(data => {
       
        const tableBody = document.querySelector('#earthquakeTable tbody');
        
        
        const listContainer = document.getElementById('earthquake-list');
        const markers = [];  





        data.forEach((data, index) => {
            const itemDiv = document.createElement('div');
            itemDiv.classList.add('item');
            itemDiv.setAttribute('data-index', index);  
    
            const numberDiv = document.createElement('div');
            numberDiv.classList.add('number');
    
            
            let color;
            if (data.magnitude < 3) {
                numberDiv.classList.add('green');
                color = 'green';
            } else if (data.magnitude < 5) {
                numberDiv.classList.add('orange');
                color = 'orange';
            } else {
                numberDiv.classList.add('red');
                color = 'red';
            }
    
            numberDiv.textContent = data.magnitude.toFixed(2);
            itemDiv.appendChild(numberDiv);
    
            const locationDiv = document.createElement('div');
            locationDiv.classList.add('location');
            locationDiv.innerHTML = `<p>${data.station}</p>`;
            itemDiv.appendChild(locationDiv);
    
            const coordinatesDiv = document.createElement('div');
            coordinatesDiv.classList.add('coordinates');
            coordinatesDiv.innerHTML = `<span class="icon">📍</span> <span> (${data.latitude}, ${data.longitude})</span>`;
            itemDiv.appendChild(coordinatesDiv);
    
            const dateDiv = document.createElement('div');
            dateDiv.classList.add('date');
            dateDiv.innerHTML = `<span class="icon">🗓️</span> <span>${data.date.substring(0,10)}</span>`;
            itemDiv.appendChild(dateDiv);

            const triggerDiv = document.createElement('div');
            triggerDiv.classList.add('date');
            triggerDiv.innerHTML = `<span class="icon">🕒</span> <span>Rel. time on: ${data.trigger_on.toFixed(2)} <br> Rel. time off: ${data.trigger_off.toFixed(2)}</span>`;
            itemDiv.appendChild(triggerDiv);
    
            listContainer.appendChild(itemDiv);
    
           
            const customIcon = {
                path: google.maps.SymbolPath.CIRCLE,
                fillColor: '#1e6091',
                fillOpacity: 1,
                scale: 10,  
                strokeWeight: 2,
                strokeColor: 'white'  
            };
    
            const marker = new google.maps.Marker({
                position: { lat: data.latitude, lng: data.longitude },
                map: map,
                icon: customIcon,  
                title: `Magnitude: ${data.magnitude}`
            });
    
            const infoWindow = new google.maps.InfoWindow({
                content: `<h4>Station: ${data.station}</h4><p>Location: ${data.latitude} ${data.longitude}</p>`,
            });
    
            markers.push({ marker, infoWindow });
    
            marker.addListener('click', () => {
                infoWindow.open(map, marker);
                map.panTo(marker.getPosition());  
            });
    
            itemDiv.addEventListener('click', () => {
                
                markers.forEach(m => m.infoWindow.close());
    
                const selectedMarker = markers[index];
                selectedMarker.infoWindow.open(map, selectedMarker.marker);
    
                map.panTo(selectedMarker.marker.getPosition());
                map.setZoom(5);
            });
        });
    })
    .catch(error => console.error('Error fetching data:', error));



}

function getNormalizedCoord(coord, zoom) {
    const y = coord.y;
    let x = coord.x;
    const tileRange = 1 << zoom;

    if (y < 0 || y >= tileRange) {
        return null;
    }

    if (x < 0 || x >= tileRange) {
        x = ((x % tileRange) + tileRange) % tileRange;
    }

    return { x: x, y: y };
}
