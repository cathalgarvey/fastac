'''nucutils - Nucleotide-wrangling functions and constants.
by Cathal Garvey
Part of the DNAmespace project. License accessible as nucutils.license.
'''
# A set of dictionaries with keys "codons", "aminos" and "starts", containing
# respectively a dictionary mapping of codons to aminos, aminos to lists of
# corresponding codons, and start codons.
from fastac import translationtables
import random

def _chunks(l, n):
    "Yield successive n-sized chunks from l."
    for i in range(0, len(l), n): yield l[i:i+n]

iupac_nucleotides = [ 'A', 'T', 'C', 'G', 'U',    # Canonical bases
                     'B', 'V', 'D', 'H',     # B=Not A, V=Not T, D=Not C, H=Not G
                     'S', 'W',               # Strong and Weak: GC vs. AT
                     'K', 'M',               # "Keto" and "aMino": GT vs. AC
                     'R', 'Y',               # "puRine" and "pYrimidine": AG vs CT
                     'N' ]         # Wildcards: any N.

aminoiupac = {'A': 'Alanine',              'B': 'Aspartic Acid or Asparagine',
              'C': 'Cysteine',             'D': 'Aspartic Acid',
              'E': 'Glutamic Acid',        'F': 'Phenylalanine',
              'G': 'Glycine',              'H': 'Histidine',
              'I': 'Isoleucine',           'K': 'Lysine',
              'L': 'Leucine',              'M': 'Methionine',
              'N': 'Asparagine',           'P': 'Proline',
              'Q': 'Glutamine',            'R': 'Arginine',
              'S': 'Serine',               'T': 'Threonine',
              'V': 'Valine',               'W': 'Tryptophan',
              'X': 'Any',                  'Y': 'Tyrosine',
              'Z': 'Glutamine or Glutamic Acid', '*': 'Stop'}

# All complements are in lowercase so that string substitution is simple:
# just set string to uppercase (if not already), replace all characters
# with lowercase complements, then if desired set resulting string to
# uppercase again.
dnacomplement = {"A":"t","T":"a","G":"c","C":"g"}
rnacomplement = {"A":"u","U":"a","G":"c","C":"g"}
dnaiupaccomplement = {   # Translates to an IUPAC sequence of potential complements.
              "A":   "t", "T": "a",
              "C":   "g", "G": "c",
              "W":   "w", "S": "s",
              "B":   "v", "V": "b",
              "H":   "d", "D": "h",
              "M":   "k", "K": "m",
              "R":   "y", "Y": "r",
              "V":   "b", "B": "v",
              "N":   "n"}
rnaiupaccomplement = {   # Translates to an IUPAC sequence of potential complements.
              "A":   "u", "U": "a", "C":   "g", "G": "c",
              "W":   "w", "S": "s", "B":   "v", "V": "b",
              "H":   "d", "D": "h", "M":   "k", "K": "m",
              "R":   "y", "Y": "r", "V":   "b", "B": "v",
              "N":   "n"}


def _uniquify(string):
    '''Reduces a string down to its component characters.
    This is a fast, order-preserving function for removing duplicates
    from a sequence/list, by "Dave Kirby"
    Found here: http://www.peterbe.com/plog/uniqifiers-benchmark'''
    seen = set()
    return [x for x in string if x not in seen and not seen.add(x)]

def deduce_alphabet(string):
    string = string.upper()
    charset = _uniquify(string)
    couldbenuc, couldbeaminos = True, True
    for character in charset:
        if character not in iupac_nucleotides:
            couldbenuc = False
        if character not in aminoiupac:
            couldbeaminos = False
    if couldbenuc:
        if "T" in charset and "U" in charset:
            if not couldbeaminos:
                raise ValueError("String contents match IUPAC code for"
                      " nucleotides but not Amino Acids, yet contains U and T"
                      " characters. Charset is '{}'".format(charset))
        elif "T" in charset:
            return "dna"
        elif "U" in charset:
            return "rna"
        else:
            # Assume DNA if ambiguous.
            return "dna"
    elif couldbeaminos:
        return "aminos"
    else:
        raise ValueError("Could not deduce an IUPAC alphabet from the input "
                         "string. Charset is '{}'".format(charset))

def get_complement_alphabet(string):
    string = string.upper()
    alphabet = deduce_alphabet(string)
    if alphabet == "aminos":
        raise ValueError("Cannot get the complement of an amino sequence.")
    elif alphabet == "rna":
        return rnaiupaccomplement
    elif alphabet == "dna":
        return dnaiupaccomplement
    else:
        raise ValueError("Could not get complement for string: "+string+"\nAlphabet deduced was: "+str(alphabet))

def get_complement(nucleotides):
    'Given a string of nucleotides (RNA *or* DNA), return reverse complement.'
    # Reverse the string:
    nucleotides = nucleotides.upper()[::-1]
    # Determine molecule type:
    basedict = get_complement_alphabet(nucleotides)
    for base in basedict.keys():
        # Replace each base with its complementary base in lowercase:
        # (lowercase means previously substituted bases aren't replaced)
        nucleotides = nucleotides.replace(base, basedict[base])
    # Return sequence in uppercase:
    return nucleotides.upper()

def translate(sequence, table, frame=1):
    sequence = sequence.upper()
    frame -= 1
    translation_table = translationtables.__dict__[table]
    aminos = []
    for codon in _chunks(sequence[frame:], 3):
        if len(codon) < 3: break
        encoded = translation_table['codons'][codon]
        aminos.append(encoded)
        if encoded == "*": break
    return ''.join(aminos)

def dumb_backtranslate(sequence, table):
    'Using the chosen table, return a back-translation of an amino sequence without codon weighting.'
    sequence = sequence.upper()
    translation_table = translationtables.__dict__[table]
    codons = []
    for amino in sequence:
        candidate_codons = translation_table['aminos'][amino]
        codons.append(random.choice(candidate_codons))
    return ''.join(codons)
