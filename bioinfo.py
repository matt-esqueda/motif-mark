 # Author: Matt Esqueda mesqueda@uoregon.edu

__version__ = "0.5"


DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')
IUPAC_bases:dict = {"A": "[Aa]", "T":"[TtUu]", "C":"[Cc]", "G":"[Gg]", "U":"[UuTt]",
              "R": "[AaGg]", "Y":"[TtCcUu]", "S":"CcGg", "W":"AaTtUu",
              "K":"[GgTtUu]", "M":"[AaCc]", "B":"[CcGgTtUu]", "D":"[AaGgTtUu]",
              "H":"[AaCcTtUu]", "V":"[AaCcGg]", "N":"[AaTtCcGgUu]"}


def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter) - 33


def qual_score(phred_score: str) -> float:
    """Calcuates the average quality score of a phred string"""
    total_score = 0
    for letter in phred_score:
        total_score += convert_phred(letter)
    return total_score/len(phred_score)


def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)


def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), 'String contains invalid characters'
    
    DNA = DNA.upper()
    Cs = DNA.count('C')
    Gs = DNA.count('G')
    return (Cs + Gs)/len(DNA)


def oneline_fasta(file):
    ''' Reads a fasta file and writes to a new fasta file with the sequence line from each record contained in one line'''
    output_file = ""
    with open(output_file, 'w') as out:                                    
        with open(file, 'r') as fa:
            for i, line in enumerate(fa):
                line = line.strip('\n')
                if i == 0:
                    print(line, file=out)
                elif line[0] == '>' and i != 0:
                    print(file=out)
                    print(line, file=out)
                else:
                    print(line, end='', file=out)
            print(file=out)


if __name__ == "__main__":

    
    assert convert_phred('!') == 0, "Convert phred does return correct phred score"
    assert convert_phred('D') == 35, "Convert phred does return correct phred score"
    print("Passed phred score tests")


    assert qual_score('ABC') == 33, "Qual score does not calculate the correct average quality score"
    assert qual_score('&&&') == 5, "Qual score does not calultate the correct average quality score"
    print("Passed quality score tests")


    assert validate_base_seq('ATCGATCG') == True, "Validate base seq does not work on DNA"
    assert validate_base_seq('GUCUUUAU', True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq('This is wrong') == False, "Validae base seq fails to recognize non-DNA"
    assert validate_base_seq('This is also wrong', True) == False, "Validate base seq fails to recognize non-RNA"
    print("Passed DNA and RNA tests")


    assert gc_content('CCGGAATT') == 0.5, "gc content does not calculate the correct GC content"
    assert gc_content('CCGGCCGG') == 1, "gc content does not calculate the correct GC content"
    assert gc_content('ATATATAT') == 0, "gc content does not calculate the correct GC content"
    print("Passes GC content tests")

    

