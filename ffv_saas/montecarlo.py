
"""
ffv.montecarlo
Calculates the fractional free volume based on random insertions

Copyright 2025 Brandon C. Tapia
"""

def distance(a, b):
    # return distance between two atoms

def run(trials):

    vdw = {"H":  1.1,
           "C":  1.7,
           "N":  1.55,
           "O":  1.5,
           "F":  1.47,
           "Si": 2.1,
           "S":  1.8,
           "Cl": 1.75,
           "Br": 1.85}

    # calculate box dimensions


    # for num_trials:
        # randomly select point within box

        # for each_atom:
            # calculate distance between point and each atom
            # if the space between the atom and the point is less than the VDW radius:
                # increment a counter
                # break

    # calculate the box volume
    # calculate the volume of points that were within the vdw radius
    # calculate the uncertainty
