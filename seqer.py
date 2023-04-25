
import doctest

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
    """ A class for proteins. Each protein has a sequence, a length, a list of all residues as words,
    and a list of all residues as abbreviations."""
    
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
        * word_seq: list<str>
        * abbrev_seq: list<str>
        
        >>> glutathione = Protein("GCE")
        "Your protein was created."
        """
        self.sequence = sequence.upper()
        self.length = len(sequence)
        self.word_seq = []
        self.abbrev_seq = []
    
    def __str__(self):
        return self.sequence.upper()
    
    def get_word_seq(self):
               
        overall_d = {"N": "Aspargine", "Q": "Glutamine", "S": "Serine", "T": "Threonine", "Y": "Tyrosine", "A": "Alanine", "G": "Glycine", "C": "Cystein", "L":"Leucine", "P": "Proline", "V": "Valine", "I": "Isoleucine", "F":"Phenylalanine", "M":"Methionine", "W": "Tryptophan", "K": "Lysine", "R": "Arginine", "H": "Histidine", "D": "Aspartic Acid", "E": "Glutamic Acid"}
        for residue in self.sequence:
            if residue in overall_d:
                self.word_seq.append(overall_d[residue])
            else:
                self.word_seq.append("("+residue + ")")
        
        
        return self.word_seq
        
    def get_abbrev_seq(self):
        overall_abbrev_d = {"N": "Asn", "Q": "Gln", "S": "Ser", "T": "Thr", "Y": "Tyr", "A": "Ala", "G": "Gly", "C": "Cys", "L":"Leu", "P": "Pro", "V": "Val", "I": "Ile", "F":"Phe", "M":"Met", "W": "Trp", "K": "Lys", "R": "Arg", "H": "His", "D": "Asp", "E": "Glu"}
        for residue in self.sequence:
            self.abbrev_seq.append(overall_abbrev_d[residue])
        
        return self.abbrev_seq
    
    def find_end_of_repetition(self, index, target_domain):
        # NEEDS TO BE FIXED SO it looks at domains and not amino accids # i think it was fixed?
        """(str, int, str) -> int
        Takes a string of amino acids, an int and a str: index and target_aa.
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
        """ Takes a seq of amino acids and returns as seq of domains with the number of amino acids in each domain.
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
    
        return domain_seq
        
        
    def has_signal_seq(self):
        """ Takes a seq of amino acids and returns true if it has a signal sequence
        >>> thing = Protein('WWWWWWWWWWWWWA')
        >>> thing.has_signalseq()
        True
        """
#         return 8 <= self.domains()[0][1] <= 16
        return 4 <= self.domains()[0][1] <= 8
    
    def has_tm_seq(self):
        """ Takes a sequence of amino acids and returns if it has a transmembrane sequence.
        >>> seq_analysis("AAK")
        """
        for domain in self.domains()[1:]: # [('Thing', 20), ('Thing 2', 20)]
#             if domain[0] == 'Non Polar' and 20 <= domain[1] <= 25:
            if domain[0] == 'Non Polar' and 10 <= domain[1] <= 12:
                return True
                
        return False
        
    def is_er_soluble(self):
        if self.has_signal_seq():
            return True
        return False
        
    def is_type_1_tm(self):
        if self.has_signal_seq() and self.has_tm_seq():
            return True
        return False
        
    def is_type_2_tm(self):
        if self.has_tm_seq():
            return True
        return False
        
        
    
#if __name__ == "__main__":
#    doctest.testmod()
    
    
    

