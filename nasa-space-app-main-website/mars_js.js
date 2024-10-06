window.onload = function() {
    initModel();
};

function initModel() {
    
    fetch('http://127.0.0.1:5000/mars-earthquake-data')
    .then(response => response.json())
    .then(data => {
        const listContainer = document.getElementById('earthquake-list');
        
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

            listContainer.appendChild(itemDiv);

            const triggerDiv = document.createElement('div');
            triggerDiv.classList.add('date');
            triggerDiv.innerHTML = `<span class="icon">🕒</span> <span>Rel. time on: ${data.trigger_on.toFixed(2)} <br> Rel. time off: ${data.trigger_off.toFixed(2)}</span>`;
            itemDiv.appendChild(triggerDiv);
    
            listContainer.appendChild(itemDiv);
        });
    })
    .catch(error => console.error('Error fetching earthquake data:', error));
}
