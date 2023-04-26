import seqer-graph as graph
import seqer-text as text
import seqer as seq

# Sample sequence: WWWWWWWWDDDDDDDDDRRRRRRRRRKK

def main():
    
    a = True
    protein_chosen = False
    
    while a:
        choice = input("Menu: (N)ew Sequence (G)raphical Analysis (T)ext Analysis (Q)uit > ")

        if choice == 'N':
                input_seq = input("Enter protein sequence: ")
                sequence = seq.Protein(input_seq)
                protein_chosen = True
                
        elif choice == 'G' and protein_chosen:
            graph.plot_charges(sequence)
            continue
        
        elif choice == 'T' and protein_chosen:
            print(text.localization(sequence))
            
            choice_2 = input("Justification? Y/N > ")
            if choice_2 == 'Y':
                text.print_justification(sequence)
                
            elif choice_2 == 'N':
                continue
            
            continue
        
        elif choice == 'Q':
            a = False
            
        else:
            print("Enter N, G, T or Q, and input a protein to analyze first.")
    
    print("Quit SEQ-ER")
