# DES_simulation
Simulation software for decay energy spectroscopy spectra given sample composition, total activity, and absorber dimensions.

## Installing this package
https://pypi.org/project/DES-simulation/0.0.1/

`pip install DES-simulation==0.0.1`

OR use `setup.py` to install  

If needed, dependencies are listed in the requirements.txt file.

`pip install -r requirements.txt`

## Running the Simulation
This software uses a command line interface with four expected arguments as shown below. The simulation then prints the results to the command window. The results include any expected escape peak values as well as the full energy peak which corresponds to the Q-value of the user-declared isotope.

### Expected Command Line Arguments
```
-i ISOTOPE, --isotope ISOTOPE
                      isotope in DES spectrum, example: Am-241
-a ACTIVITY, --activity ACTIVITY
                      activity of sample in Bq
-x WIDTH, --width WIDTH
                      absorber width in cm
-y LENGTH, --length LENGTH
                      absorber length in cm
-z THICKNESS, --thickness THICKNESS
                      absorber thickness in cm
```

## Brief example
Command line input:

`python DES_simulation.py -i Am-241 -a 1.0 -x 0.5 -y 0.5 -z 0.10`

The output for this command line input is:
```
Am-241 will have escape peaks at:  [5623.9212 5578.2803]
Am-241 has a full energy peak at:  5637.8212
```

**Output energy is always in units of ***keV!*****

See Expected Command Line Arguments Section above for information on the expected units for user input values.
