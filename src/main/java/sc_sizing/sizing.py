import json
import numpy as np

## Constants
G = 6.674e-11

M_sun = 1.989e30
M_earth = 5.972e24
M_mars = 6.39e23
M_saturn = 5.6834e26
M_titan = 1345.5e20
M_moon = 7.34767309e22

mu_sun = G*M_sun
#mu_earth = G*M_earth
mu_earth = 3.986004415e14
mu_mars = G*M_mars
mu_saturn  = G*M_saturn
mu_titan = G*M_titan
mu_moon = G*M_moon

R_earth = 6378.1363
R_sun = 695508
R_mars = 3389
R_saturn = 58232
R_titan= 2575
R_moon = 1737

AU = 1.496e11
a_earth = 1*AU
a_mars = 1.524*AU
a_saturn= 9.557*AU
a_moon = 0.3844e9

#J2_earth = 1.08262668e-3
J2_earth = 1.0826362e-3

def get_instrument_lists(file_name):
    # -Returns the Instrument lists from input JSON-
    # Open file
    filePath = "./inputs/" + file_name
    with open('./inputs/test_input.json') as f:
        input_data = json.load(f)

    # Read every satellite in space segment
    instrument_lists = []
    for i in range(len(input_data['spaceSegment'][0]['satellites'])):
        tempList = []
        for j in range(len(input_data['spaceSegment'][0]['satellites'][i]['payload'])):
            tempList.append(input_data['spaceSegment'][0]['satellites'][i]['payload'][j]['acronym'])
        instrument_lists.append(tempList)

    return instrument_lists

def get_orbit_lists(file_name):
    # -Returns the Orbit lists from input JSON-
    # Open file
    filePath = "./inputs/" + file_name
    with open('./inputs/test_input.json') as f:
        input_data = json.load(f)

    # Read every satellite in space segment
    orbit_lists = []
    for i in range(len(input_data['spaceSegment'][0]['satellites'])):
        tempList = []
        # Translate inputs into VASSAR format
        orbit_data = input_data['spaceSegment'][0]['satellites'][i]['orbit']
        tempList.append( translate_orbit(orbit_data) )
        orbit_lists.append(tempList)

    return orbit_lists

def translate_orbit(orbit_data):
    # -Translate inputs into VASSAR format-
    # Initialize result and unpackage input
    type = ""
    h = 0.0
    inc = ""

    a = orbit_data['semimajorAxis']
    i = orbit_data['inclination']
    e = orbit_data['eccentricity']
    arg = orbit_data['periapsisArgument']
    raan = orbit_data['rightAscensionAscendingNode']
    anom = orbit_data['trueAnomaly']
    epoch = orbit_data['epoch']
    time = orbit_data['time']

    if e <= 0.1:
        h = int( a - R_earth )
    else:
        h = int( (a - a * e) - R_earth )

    if is_sso(a,e,i):
        type = "SSO"
        inc = "SSO"
    elif is_geo(a):
        type = "GEO"
        if abs(i - 90) <= 0.1:
            inc = "polar"
        else:
            inc = "NA"
    else:
        type = "LEO"
        time = "NA"
        if abs(i - 90) <= 0.1:
            inc = "polar"
        else:
            inc = "NA"

    return type + "-" + str(h) + "-" + inc + "-" + time

def is_sso(a,e,i):
    n = np.sqrt( mu_earth / np.power(a,3) )
    p = a*(1 - np.power(e,2))
    i_sso = (180/np.pi) * np.arccos( -1.227*10e-4*(1/n)*( np.power(p,2) / np.power(R_earth, 2)) )

    return abs(i - i_sso) <= 0.01

def is_geo(a):
    T = 2*np.pi*np.sqrt( np.power(a,3) / mu_earth )
    return  abs(T - 24*3600) <= 60