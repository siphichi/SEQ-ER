
import doctest
import matplotlib as mp
import numpy as np

"""
polar_d = {"N": "Aspargine", "Q": "Glutamine", "S": "Serine", "T": "Threonine", "Y": "Tyrosine"}
        nonpolar_d = {"A": "Alanine", "G": "Glycine", "C": "Cystein", "L":"Leucine", "P": "Proline", "V": "Valine", "I": "Isoleucine", "F":"Phenylalanine", "M":"Methionine", "W": "Tryptophan"}
        basic_d = {"K": "Lysine", "R": "Arginine", "H": "Histidine"}
        acidic_d = {"D": "Aspartic Acid", "E": "Glutamic Acid"}
"""

def find(aa):
    """ Takes an amino acid str and finds its domain.
    """
    domain = ""
    polar = "NQSTYG"
    nonpolar = "AGCLPVIFMW"
    basic = "KRHG"
    acidic = "DEG"
    
    if aa in polar:
        domain = "Polar"
    elif aa in nonpolar:
        domain = "Non Polar"
    elif aa in basic:
        domain = "Basic"
    elif aa in acidic:
        domain = "Acidic"
    
    return domain
    

class Protein:
    """ A class for proteins. Each protein has a sequence, a length, and a list of domains."""
    
    polar = "NQSTY"
    nonpolar = "AGCLPVIFMW"
    basic = "KRH"
    acidic = "DE"

    def __init__(self, sequence):
        """(str) -> Protein
        A constructor that takes a protein sequence string as explicit input and returns nothing.
        Creates the following instance attributes:
        ∗ sequence: str
        ∗ length: int
        * domains: list<str>
        
        >>> glutathione = Protein("GCE")
        """
        self.sequence = sequence.upper()
        self.length = len(sequence)
        self.domain_list = []
    
    def __str__(self):
        return self.sequence.upper()
    
    
    def find_end_of_repetition(self, index, target_domain):
        """(str, int, str) -> int
        Takes a string of amino acids, an int index and a target sequence target_aa.
        Returns the index of the last consecutive occurrence of the target_aa. 
        >>> find_end_of_repetition("AAAAADDEE", 2, "A")
        4
        >>> find_end_of_repetition([1, 2, 3, 4, 5, 6, 7], 7, 7)
        6
        >>> find_end_of_repetition([1, 1, 1, 1, 1], 0, 1)
        4
        >>> find_end_of_repetition([1], 0, 1)
        0
        """
            
        while index < len(self.sequence):
            if find(self.sequence[index]) == target_domain:
                index += 1
            else:
                break
                
        return index - 1
    
    def domains(self):
        """ Takes a sequence of amino acids and returns as seq of domains with the number of amino acids in each domain.
        The domain can be polar, acidic, basic, and hydrophobic.
        """
        domain_seq = []
        i = 0
        previous_end_of_rep=0
    

        while i < len(self.sequence):
            i = self.find_end_of_repetition(i, find(self.sequence[i]))
            domain_seq += [ ( find(str(self.sequence[i])), int(i + 1 - previous_end_of_rep) ) ]
            i += 1
            previous_end_of_rep = i
    
        self.domain_list = domain_seq
        return domain_seq
        
        
    def has_signal_seq(self):
        """ Takes a seq of amino acids and returns true if it has a signal sequence
        >>> protein = Protein('WWWWWWWWWWWWWA')
        >>> protein.has_signal_seq()
        True
        """

        return 8 <= self.domains()[0][1] and self.domains()[0][0] == 'Non Polar'
    
    
    def has_tm_seq(self):
        """ Takes a sequence of amino acids and returns if it has a transmembrane sequence.
        >>> protein = Protein('WWWWWWWWWWWWWA')
        >>> protein.has_tm_seq()
        False
        """
        for domain in self.domains()[1:]: # [('Domain', 20), ('Domain 2', 10)]
            if domain[0] == 'Non Polar' and 20 <= domain[1]:
                return True
                
        return False
        
    def has_er_motif(self):
        """ Takes a sequence of amino acids and returns if it has the KDEL sequence
            at its C terminus."""
        
        return "KDEL" in self.sequence[-4:]
    
    def has_tm1_er_motif(self):
        """ Takes a sequence of amino acids and returns if it has the KK** sequence
            at its C terminus."""
        
        return "KK" in self.sequence[-2:]
        
    def has_tm2_er_motif(self):
        """ Takes a sequence of amino acids and returns if it has the M**RR sequence
            at its N terminus."""
        
        return self.sequence[0] == "M" and self.sequence[3:5] == "RR"
        
    def has_golgi_motif(self):
        """ Takes a sequence of amino acids and returns if it has the FF sequence
            at its C terminus."""
        
        return "FF" in self.sequence[-2:]
        
    def net_charge(self):
        """ Takes a sequence of amino acids and calculates their overall charge."""

        charge = 0.0
        
        for i in range(len(self.sequence)):
            pk_list = {"K": 10.5, "R": 12.5, "H": 6, "D": 3.9, "E": 4.3, "S": 13, "T": 13, "C": 8.3}
            
            if self.sequence[i] in pk_list:
                pk = pk_list[self.sequence[i]]
                
                if self.sequence[i] in "KRH":
                    charge += (1 - 1/(10**(pk-7.2)+1))
                    
                elif self.sequence[i] in "DESTC":
                    charge = charge - (1/(10**(pk-7.2)+1))
                    
        return charge/(i+1)
    
        # ER = more oxidizing
    def partial_charge(self, start_i = 0, end_i = -1):
        """Takes a sequence of amino acids and calculates the partial charge of a domain, only in terms of charged residues."""
        
        charge = 0.0
        n=0
        
        for i in range(len(self.sequence[start_i : end_i])+1):
            pk_list = {"K": 10.5, "R": 12.5, "H": 6, "D": 3.9, "E": 4.3, "S": 13, "T": 13, "C": 8.3}

            if self.sequence[i] in pk_list:
                
                n +=1
                pk = pk_list[self.sequence[i]]
                
                if self.sequence[i] in "KRH":
                    charge += (1 - 1/(10**(pk-7.2)+1))
                    
                elif self.sequence[i] in "DESTC":
                    charge = charge - (1/(10**(pk-7.2)+1))
                    
        print(charge)
        
        if n>0:
            partial_charge = charge/n
        else:
            partial_charge = 0
            
        return partial_charge
        
        
    def indiv_charge(self, index):
        """Takes a sequence of amino acids and calculates the partial charge of a domain, only in terms of charged residues."""
        
        charge = 0.0
        
        pk_list = {"K": 10.5, "R": 12.5, "H": 6, "D": 3.9, "E": 4.3, "S": 13, "T": 13, "C": 8.3}

        if self.sequence[index] in pk_list:
            pk = pk_list[self.sequence[index]]
            
            if self.sequence[index] in "KRH":
                charge += (1 - 1/(10**(pk-7.2)+1))
                
            elif self.sequence[index] in "DESTC":
                    charge = charge - (1/(10**(pk-7.2)+1))
                    
        return charge
        
        
        
        
#if __name__ == "__main__":
#    doctest.testmod()
