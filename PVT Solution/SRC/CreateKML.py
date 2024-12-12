#!/usr/bin/env python

########################################################################
# CreateKML.py:
# It produces a KML file with the Long, Lat, and Alt of LEO RCVR
#
#  Project:        SENTUS
#  File:           CreateKML.py
#
#   Author: Nerea Sánchez
#   Copyright 2024 GNSS Academy / Nerea Sánchez
#
# Usage:
#   CreateKML.py $SCEN_PATH
########################################################################

import simplekml
import sys
import os
from InputOutput import PosIdx 

def createKML(file_path, Scen):

    '''
    Create a KML file with coordinates from a .dat file.

    Parameters
    ----------
    file_path : str
        The path to the input .dat file.
    Scen : str
        The base path where the output KML file will be saved.
    '''

    # Check if the .dat file exists
    if not os.path.exists(file_path):
        print(f"Error: the file {file_path} does not exist.")
        sys.exit(1)

    # Step 1: Read the .dat file and extract the columns
    coordinates = []  # List to store the coordinates (longitude, latitude, altitude)

    # Read the .dat file
    with open(file_path, 'r') as file:
        for line in file:
            # Skip lines that don't contain data
            if not line.strip():
                continue
            
            # Split the values in each line (assuming they are space-separated)
            values = line.split()

            try:
                # Extract the columns using PosIdx for longitude, latitude, and altitude
                lon = float(values[PosIdx["LONEST"]])  # Longitude
                lat = float(values[PosIdx["LATEST"]])  # Latitude
                alt = float(values[PosIdx["ALTEST"]])  # Altitude

                # Add the coordinates (longitude, latitude, altitude) to the list
                coordinates.append((lon, lat, alt))
            except (IndexError, ValueError) as e:
                continue

    # Step 2: Create the KML file with simplekml
    kml = simplekml.Kml()

    # Add a LineString to represent the satellite's orbit track
    linestring = kml.newlinestring(name="Satellite Orbit Track")

    # Add coordinates (longitude, latitude, altitude) to the LineString
    linestring.coords = coordinates

    # Set the altitude mode to absolute so the line represents the orbit track in space
    linestring.altitudemode = simplekml.AltitudeMode.absolute  # Altitude relative to sea level
    linestring.extrude = 1  # Optional: This makes the line extrude towards the Earth's surface

    # Step 3: Save the KML file
    output_file = os.path.join(Scen, 'OUT', 'PVTS', 'LEOCoord.kml')
    kml.save(output_file)

    print(f"The KML file '{output_file}' has been successfully created.")
