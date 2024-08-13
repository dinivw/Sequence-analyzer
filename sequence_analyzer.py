'''
Final project Assignment - Sequence analyzer
Name - Dinithi Wanniarachchi
Index_no - s14720
'''

#importing the necessary packages
import numpy as np
import matplotlib.pyplot as plt

#Main class Sequence
class Sequence:
    def __init__(self,sequence):
        self.sequence = sequence
        self.GCcontent = self.getGCcontent()


    #static method to get only the needed sequence out of a fasta file
    #this is an additional method which is not mentioned in the assignment
    @staticmethod
    def fastasplit(file):

        #empty string variable to save the sequence
        sequence = ''

        #open the fasta file and remove unwanted characters
        with open(file) as file:
            for line in file:
                if line != '\n':
                    line = line.strip()
                    #if '>' not in line add each line to the empty variable created
                    if ">" not in line:
                        sequence += line

        return sequence

    #a static method to get the type of any sequence (DNA,RNA,protein) that is inputed in to the program
    @staticmethod
    def getSeqType(sequence):

        #save the protein except for A,C,G,T in a list
        proteinlist = ["M", "N", "R", "S", "I", "Q", "H", "P", "R", "L", "E", "D", "V", "O", "Y", "W", "L", "F", "W"]
        # for each protein in the list if the protein is in sequence it is a protein sequence
        # if U present in the sequence it is a mRNA sequence
        # else U is absent in the sequence take it as a DNA sequence
        for protein in proteinlist:
            if protein in sequence:
                return 'Protein'
                break

            elif "U" in sequence:
                return 'RNA'
                break

            else:
                return 'DNA'
                break

    #an instance method to get the GC content of a sequence
    def getGCcontent(self):

        #set a variable equal to 0 to store the G and C count
        GCcount = 0
        #for each base in sequence if the base is G or C increase the count by one
        for base in self.sequence:
            if base == 'G' or base == 'C':
                GCcount += 1

        #get the GC content by dividing the count by length of the sequence
        GCcontent = GCcount / len(self.sequence)
        return GCcontent

    #an instance method to get all possible combinations of K-mers when k value is given,for a sequence

    def getKmers(self,K):

        #a dictionary to save all possible combinations and their frequency
        counts = {}

        #number of kmers are found by substracting K value from length of the sequence and adding one
        #this makes sure to iterate over all possible Kmer substrings
        num_kmers = len(self.sequence) - K + 1

        for i in range(num_kmers):
            kmer = self.sequence[i:i + K]
            # Add the kmer to the dictionary if it's not present
            if kmer not in counts:
                counts[kmer] = 0
            #increase the value accroding to the frequency of specific Kmer being found
            counts[kmer] += 1

        return counts

#subclass for DNA sequences
class DNAsequence(Sequence):
    #constructor method
    def __init__(self,sequence):
        super().__init__(sequence)
        self.sequence = sequence
        self.Basecount = self.getBaseCount()
        self.ATcontent = self.getATcontent()
        self.GCcontent = self.getGCcontent()
        self.Transcribed_seq = self.getTranscribedseq()

    #An instance method to get th count of each base in a neucleotide sequence
    #this method returns the count of each base and also a histogram comparing base counts
    def getBaseCount(self):

        #set variable for each base equal to 0
        adenineCount = 0
        ThymineCount = 0
        GuanineCount = 0
        CytosineCount = 0

        #for each base in the sequence increase the counts if the respective base is found
        for base in self.sequence:
            if base == 'A':
                adenineCount += 1
            if base == 'T':
                ThymineCount += 1
            if base == 'G':
                GuanineCount += 1
            if base == 'C':
                CytosineCount += 1

        #A histrogram to comapre the counts of each base

        #each count name and the value is stored in a dictionary
        data = {'Adenine Count': adenineCount, 'Thymine Count': ThymineCount, 'Guanine count': GuanineCount,'Cytosine Count': CytosineCount}

        #the names and counts are extracted as Bases and values
        Base = list(data.keys())
        value = list(data.values())

        #then the bases and their values are plotted in the graph
        fig = plt.figure(figsize=(10, 5))
        plt.bar(Base, value, color='maroon', width=0.4)
        plt.xlabel('Base')
        plt.ylabel('Base count')
        plt.title('Count of each base in the DNA sequence')
        plt.savefig('Base counts - DNA.jpg')

        return ('Adenine count:',adenineCount ,'Thymine count:',ThymineCount,'Guanine count:',GuanineCount,'Cytosine count:',CytosineCount)

    #an instance method to get the AT content of a sequence
    def getATcontent(self):

        #set a variable equal to 0 to save the count
        ATcount = 0
        #for each base in the sequence if the base is A or T increase the count by 1
        for base in self.sequence:
            if base == 'A' or base == 'T':
                ATcount += 1

        #get the AT content by diving the count from length of the sequence
        ATcontent = ATcount / len(self.sequence)
        return ATcontent


    #a class method to compare AT and GC counts of two DNA sequences
    #this method outputs a graphical image with comparisons
    @classmethod
    def compareATandGCcontents(cls,sequence1, sequence2):
        #define emty variables for the counts and set them equal to zero
        ATcount1 = 0
        ATcount2 = 0
        GCcount1 = 0
        GCcount2 = 0

        #GO through each base in the sequence find the AT and GC counts in both the sequences
        for base in sequence1:
            if base == 'A' or base == 'T':
                ATcount1 += 1
            if base == 'G' or base == 'C':
                GCcount1 += 1

        for base in sequence2:
            if base == 'A' or base == 'T':
                ATcount2 += 1
            if base == 'G' or base == 'C':
                GCcount2 += 1

        #calculate the AT and GC contents by dividing the counts from length of the sequence
        ATcontent1 = ATcount1 / len(sequence1)
        ATcontent2 = ATcount2 / len(sequence2)
        GCcontent1 = GCcount1 / len(sequence1)
        GCcontent2 = GCcount2 / len(sequence2)

        #save the AT and GC contents of both sequences in lists
        seq1 = [ATcontent1, GCcontent1]
        seq2 = [ATcontent2, GCcontent2]

        #get an array of 0 and 1
        r = np.arange(2)
        #define the width of the bar
        width = 0.15

        #plot the bars graph
        plt.bar(r, seq1, color='b', width=width, edgecolor='black', label='Sequence 1')
        plt.bar(r + width, seq2, color='g', width=width, edgecolor='black', label='Sequence 2')

        #label axis and titles
        plt.xlabel('Content type')
        plt.ylabel('Content value')
        plt.title('AT and GC content of sequences')

        #get the labels in proper manner
        plt.xticks(r + width / 2, ['AT content', 'GC content'])
        plt.legend()

        #save the figure
        plt.savefig('Comparison graph.jpg')

    #an instance method to transcribe the DNA sequence and get the RNA sequence
    def getTranscribedseq(self):

        #an empty string variable to save the sequence
        transcribed_seq = ''

        #for each base in sequence if the base is T it is converted to U, or else left as it is
        for base in self.sequence:
            if base == 'T':
                transcribed_seq += 'U'
            else:
                transcribed_seq += base

        return transcribed_seq

#a subclass for RNA sequences
class RNAseq(Sequence):
    #constructor method
    def __init__(self,sequence):
        super().__init__(sequence)
        self.sequence = sequence
        self.Basecount = self.getBaseCounts()
        self.ATcontent = self.getATcontent()
        self.GCcontent = self.getGCcontent()
        self.TranslatedSeq = self.getTranslatedSeq()

    #a method to get the base counts of a RNA sequence
    #this method outputs the counts of each base and a graph comapring the base counts
    def getBaseCounts(self):

        #define new variables to sabe the base counts and set them equal to zero
        adenineCount = 0
        UracilCount = 0
        GuanineCount = 0
        CytosineCount = 0

        #go through each base in the sequence and get the base counts
        for base in self.sequence:
            if base == 'A':
                adenineCount += 1
            if base == 'U':
                UracilCount += 1
            if base == 'G':
                GuanineCount += 1
            if base == 'C':
                CytosineCount += 1

        # A histrogram to comapre the counts of each base

        # each count name and the value is stored in a dictionary
        data = {'Adenine Count': adenineCount, 'Uracil Count': UracilCount, 'Guanine count': GuanineCount,
                'Cytosine Count': CytosineCount}

        # the names and counts are extracted as Bases and values
        Base = list(data.keys())
        value = list(data.values())

        # then the bases and their values are plotted in the graph
        fig = plt.figure(figsize=(10, 5))
        plt.bar(Base, value, color='maroon', width=0.4)
        plt.xlabel('Base')
        plt.ylabel('Base count')
        plt.title('Count of each base in the sequence')
        plt.savefig('Base counts-RNA.jpg')

        return ('Adenine count:', adenineCount,'Uracil count:', UracilCount,'Guanine count:', GuanineCount,'Cytosine count:', CytosineCount)

    #a method to get the AT content of a RNA sequence
    def getATcontent(self):
        #set a new variable equal to 0 to store the AT content
        AUcount = 0
        #for each base in the sequence if the base is A ou U increase the count by one
        for base in self.sequence:
            if base == 'A' or base == 'U':
                AUcount += 1

        #get the AU content by dividing the count from length of the sequence
        ATcontent = AUcount / len(self.sequence)
        return ATcontent

        # a class method to compare AT and GC counts of two DNA sequences
        # this method outputs a graphical image with comparisons

    # a class method to compare AT and GC counts of two DNA sequences
    # this method outputs a graphical image with comparisons
    @classmethod
    def compareATandGCcontents(cls, sequence1, sequence2):
        # define emty variables for the counts and set them equal to zero
        ATcount1 = 0
        ATcount2 = 0
        GCcount1 = 0
        GCcount2 = 0

        # Go through each base in the sequence find the AT and GC counts in both the sequences
        for base in sequence1:
            if base == 'A' or base == 'U':
                ATcount1 += 1
            if base == 'G' or base == 'C':
                GCcount1 += 1

        for base in sequence2:
            if base == 'A' or base == 'U':
                ATcount2 += 1
            if base == 'G' or base == 'C':
                GCcount2 += 1

        # calculate the AT and GC contents by dividing the counts from length of the sequence
        ATcontent1 = ATcount1 / len(sequence1)
        ATcontent2 = ATcount2 / len(sequence2)
        GCcontent1 = GCcount1 / len(sequence1)
        GCcontent2 = GCcount2 / len(sequence2)

        # save the AT and GC contents of both sequences in lists
        seq1 = [ATcontent1, GCcontent1]
        seq2 = [ATcontent2, GCcontent2]

        # get an array of 0 and 1
        r = np.arange(2)
        # define the width of the bar
        width = 0.15

        #plot the bar graph
        plt.bar(r, seq1, color='b', width=width, edgecolor='black', label='Sequence 1')
        plt.bar(r + width, seq2, color='g', width=width, edgecolor='black', label='Sequence 2')

        # label axis and titles
        plt.xlabel('Content type')
        plt.ylabel('Content value')
        plt.title('AT and GC content of sequences')

        # get the labels in proper manner
        plt.xticks(r + width / 2, ['AU content', 'GC content'])
        plt.legend()

        #save the graph
        plt.savefig('Comparison graph.jpg')


    #an instance method to get the translated amino acid sequence from a RNA sequence
    def getTranslatedSeq(self):
        #define an empty string variable to save the output sequence
        translated_seq = ''

        #store the codons and their respective amino acids in a codon dictionary
        codonDic = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T', 'AGA': 'R',
         'AGC': 'S', 'AGG': 'R', 'AGU': 'S', 'AUA': 'I', 'AUC': 'I', 'AUG': 'M', 'AUU': 'I', 'CAA': 'Q', 'CAC': 'H',
         'CAG': 'Q', 'CAU': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R',
         'CGU': 'R', 'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAU': 'D',
         'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G', 'GUA': 'V',
         'GUC': 'V', 'GUG': 'V', 'GUU': 'V', 'UAA': 'O', 'UAC': 'Y', 'UAG': 'O', 'UAU': 'Y', 'UCA': 'S', 'UCC': 'S',
         'UCG': 'S', 'UCU': 'S', 'UGA': 'O', 'UGC': 'C', 'UGG': 'W', 'UGU': 'C', 'UUA': 'L', 'UUC': 'F', 'UUG': 'L',
         'UUU': 'F'}

        #go through the sequence and replace the codons with their respective amino acids using the dictionary
        for i in range(0,len(self.sequence)-3,3):
            codon = self.sequence[i:i + 3]

            #stop the translation when a stop codon is met
            if codon == 'UAA' or codon == 'UAG' or codon == 'UGA':
                break
            translated_seq += codonDic[codon]

        return translated_seq

#subclass for protein sequences
class ProteinSeq(Sequence):
    #constructor method
    def __init__(self,sequence):
        super().__init__(sequence)
        self.sequence = sequence
        self.Hydrophobicity = self.getHydrophobicity()
        self.MolecularWeight = self.getMolecularWeight()

    #an instance method to get the hydrophobicity of a specific sequence
    def getHydrophobicity(self):

        #define an empty variable to store the count of hydrophobic amino acids and set it equal to zero
        count = 0

        #a list with the hydrophobic amino acids stored in it
        aminoacidlist = ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V']

        #for each protein in the list and for each base in sequence, if the protein is same as the base increase the count by one
        for protein in aminoacidlist:
            for base in self.sequence:
                if base == protein:
                    count += 1

        #get the hydrophobicity by dividing the count from length of the sequence
        hydrphobicity = count / len(self.sequence)
        return hydrphobicity

    #a method to get the molecular weight of an amino acid sequence
    def getMolecularWeight(self):

        #a dictionary which contains all the amino acids and their specific molecular weight
        protein_weights = {'A': 89.0932, 'C': 121.1582, 'D': 133.1027, 'E': 147.1293, 'F': 165.1891, 'G': 75.066,
                           'H': 155.1546, 'I': 131.1729, 'K': 146.1876, 'L': 131.1729, 'M': 149.2113, 'N': 132.1179,
                           'O': 255.3134, 'P': 115.1305, 'Q': 146.1445, 'R': 174.2017, 'S': 105.0926, 'T': 119.1192,
                           'U': 168.0532, 'V': 117.1463, 'W': 204.2252, 'Y': 181.1885}

        #define a variable to store the molecular weight and set it equal to zero
        molecular_weight = 0

        #for each amino acid in the sequence get the specific molecular weight from the dictionary and obtain total weight
        for aminoacid in self.sequence:
            molecular_weight += protein_weights[aminoacid]
            #reduce the molecular weight of water since when peptide bonds are formed in between the amino acid a water molecule us removed
        total_molecular_weight_water = (len(self.sequence) - 1) * 18.01
        molecular_w = molecular_weight - total_molecular_weight_water

        return molecular_w


