from Bio import SeqIO

# Load the DNA seqence from "sequence.txt" using BioPython
fname = "sequence.txt"
dna = str(SeqIO.read(fname, "fasta").seq)
#print(dna) dna string is the sense strand
dnacomp = str(dna)
dnacomp = dnacomp.replace('a','x')
dnacomp = dnacomp.replace('t','a')
dnacomp = dnacomp.replace('x','t')
dnacomp = dnacomp.replace('g','y')
dnacomp = dnacomp.replace('c','g')
dnacomp = dnacomp.replace('y','c')
dnacomp = dnacomp[::-1]
#print(dnacomp) dnacomp string is the antisense strand

# Initialise a dictionary to store codons
codons = {}

# Read in the codons and amino acids from "codons.txt" and store in
# the dictionary.
in_file = open('codons.txt')
contents = in_file.read()
in_file.close()
contents = contents.split()
#assign first and then every other string (codon) as the key
#i is codon, i+1 is the value (amino acid)
for i in range(0,len(contents),2):
    codons[contents[i]]=contents[i+1]

# Initialise your translation
orf1 = ""

for count in range(0,3):
    orf1 += "\nORF %s \n" % str(count+1)
    #%s and str(count+1) is so at the beginning of each translation it prints
    #'ORF ' then the number 
    for i in range(count, len(dna), 3):
        #count 0, 1, 2, starts read at different nucleotides, so new reading frames
        codon = dna[i:i+3]
        if len(codon) == 3:
            orf1 += codons[codon]
for count in range(0,3):
    orf1 += "\nORF %s \n" % str(count+4)
    for i in range(count, len(dnacomp), 3):
        codon = dnacomp[i:i+3]
        if len(codon) == 3:
            orf1 += codons[codon]
print(orf1)
