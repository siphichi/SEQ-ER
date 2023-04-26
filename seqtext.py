import numpy as np
from seqer import *


def analyze_signals(sequence):
    """ Analyzes a protein sequence based on signal sequences and motifs."""

    
    tm1 = 0
    tm2 = 0
    er_soluble = 0
    cytosol_soluble = 1
    golgi = 0
    
    if sequence.has_signal_seq():
        tm1 += 1
        tm2 += 1
        er_soluble += 1
        cytosol_soluble -=1
        golgi += 0.2
        
    if sequence.has_tm_seq():
        tm1 += 1
        tm2 += 1
        er_soluble -= 1
        cytosol_soluble -=1
        golgi += 0.2
        
    if sequence.has_er_motif():
        tm1 += 1
        tm2 += 1
        er_soluble += 1
        golgi += 0.2
    
    if sequence.has_tm1_er_motif():
        tm1 += 1
        tm2 -= 1
        er_soluble -= 1
        golgi += 0.2
    
    if sequence.has_tm2_er_motif():
        tm1 -= 1
        tm2 += 1
        er_soluble -= 1
        golgi += 0.2
            
    if sequence.has_golgi_motif():
        golgi = 0.2
    
    results = [tm1/4, tm2/4, er_soluble/2, cytosol_soluble, golgi/1]
    return results
    
        
def analyze_charge(sequence):
    """ Analyzes the net charge and domain charges of a protein sequence"""
    
    # Could be adjusted to take into account placement of TM areas
    
    net = sequence.net_charge()
    polar = 0
    non_polar = 0
    basic = 0
    acidic = 0
    
    for i in range(1, len(sequence.domains())):
        if sequence.domains()[i][0] == 'Polar' and sequence.domains()[i][1] > 5:
            polar += 1
        elif sequence.domains()[i][0] == 'Non Polar' and sequence.domains()[i][1] > 5:
            non_polar += 1
        elif sequence.domains()[i][0] == 'Basic' and sequence.domains()[i][1] > 5:
            basic += 1
        elif sequence.domains()[i][0] == 'Acidic' and sequence.domains()[i][1] > 5:
            acidic += 1
            
    polarity = polar/i
    hydrophobicity = non_polar/i
    basicity = basic/i
    acidity = acidic/i
    
    results = [polarity, hydrophobicity, basicity, acidity, net]
    return results
    
def localization(sequence):
    signals = analyze_signals(sequence)
    charges = analyze_charge(sequence)
    likely_localization = "Unknown"
    
    # TM1
    if signals.index(max(signals)) == 0: # and charges[1] > 2 :
        likely_localization = "ER Type 1 Transmembrane"
        
    # TM2
    elif signals.index(max(signals)) == 1:# and charges[1] > 2:
        likely_localization = "ER Type 2 Transmembrane"
        
    # ER
    elif signals.index(max(signals)) == 2 : #and charges[3] > charges[2] and charges[4] < 0:
        likely_localization = "ER Soluble"
        
    # Cytosol
    elif signals.index(max(signals)) == 3 : #and charges[2] > charges[3] and charges[4] > 0:
        likely_localization = "Cytosol Soluble"
    
    # Golgi
    elif signals.index(max(signals)) == 4  :# and charges[1] > 2:
        likely_localization = "Golgi Soluble"

    return likely_localization

def print_justification(sequence):
    signals = analyze_signals(sequence)
    charges = analyze_charge(sequence)
    
    print("Chance of TM 1: " + str(signals[0]*100)+ "%")
    print("Chance of TM 2: " + str(signals[1]*100)+ "%")
    print("Chance of ER Soluble: " + str(signals[2]*100)+ "%")
    print("Chance of Cytosol Soluble: " + str(signals[3]*100) + "%")
    print("Change of Golgi: "+ str(signals[4]*100) + "%")
    print("--------------------------------------")
    print("Polarity: " + str(charges[0]*100) + "%")
    print("Hydrophobicity: " + str(charges[1]*100) + "%")
    print("Basicity: " + str(charges[2]*100) + "%")
    print("Acidity: " + str(charges[3]*100) + "%")
    print("Net charge: " + str(charges[4]))
    print("--------------------------------------")
    print(sequence.domains())
    
    
    
    
    