Our challenge is to develop a computer program capable of analyzing real seismic data from the Apollo missions and the Mars InSight Lander to identify seismic quakes hidden within the noisy recordings. In planetary seismology, missions often struggle with the power requirements needed to send continuous, high-resolution seismic data back to Earth. Since transmitting large amounts of data over vast distances requires significant energy, only a fraction of the collected data is typically scientifically valuable. This makes it crucial to devise a method that can distinguish between useful signals and noise, allowing only important seismic data to be transmitted.

Our initial approach focused on using machine learning models and pattern recognition techniques. However, we quickly realized that there was not enough data to support this strategy due to the differences in the seismic characteristics between Earth, the Moon, and Mars. With varying frequencies and velocities between these planetary bodies, a direct application of AI models proved unfeasible. Thus, we pivoted towards a more traditional approach involving signal processing.

Our solution involves the following steps:

    Data Loading and Filtering: We begin by loading the raw seismic data and applying a bandpass filter to isolate the frequencies of interest.
    Signal Energy Analysis: To get a clearer view of potential events, we calculate the Root Mean Square (RMS) of the filtered data, which gives us an understanding of the energy distribution.
    Event Detection Using Triggers: We define trigger_on and trigger_off thresholds to detect the start and end of seismic events. Initially, these thresholds are set to static values, but we also incorporate adaptive dynamic triggers if the initial settings do not effectively capture the quakes.
    Impulse Noise Mitigation: To avoid picking up unwanted spikes or false positives, we implemented logic to filter out impulse spikes that slip through the bandpass filter.
    Visualization: Finally, we plot various graphs to visually confirm the detection process:
        Raw seismic data
        Filtered seismic data
        RMS energy with defined thresholds
        Raw data overlaid with detected triggers.

This approach has allowed us to efficiently navigate the complexities of planetary seismic data, where low signal-to-noise ratios and unique noise patterns present significant challenges. While our strategy relies on signal processing, we’re continually refining our method to improve its robustness and adaptability to new datasets.

