'''
Final project assignment
Name-Dinithi Wanniarachchi
Index_no - s14720
Main program
'''

#import the python file containg all the classes and methods
import sequence_analyzer

print('-------------------------------------------Sequence Analyzer---------------------------------------------\n')
#this program takes a maximum of two files only
#the files are saved in a list
files = []
while len(files) <1 or len(files)>2:
    #input the needed files
    files = input('Input one or two fasta files separated by comma:\n').split(',')

#if only one file is inputed relevant methods can be done using the sequence of that file
if len(files) == 1:
    file_1 = files[0]
    print('\nYou entered one file:',file_1,'\n')
    #extracting the sequence from the fasta file
    sequence = sequence_analyzer.Sequence.fastasplit(file_1)

    #file type is obtained
    sequence_type_file1 = sequence_analyzer.Sequence.getSeqType(sequence)
    # print(sequence_type_file1)

    #based on the sequence type an instance of that sequence is made
    if sequence_type_file1 == 'DNA':
        object = sequence_analyzer.DNAsequence(sequence)

    if sequence_type_file1 == 'RNA':
        object = sequence_analyzer.RNAseq(sequence)

    if sequence_type_file1 == 'Protein':
        object = sequence_analyzer.ProteinSeq(sequence)

    #User can input the operation that needs to be extracted from the sequence out of the listed methods
    Operation = input(
        "Type the operation:\n" "(DNA:Seqtype,Kmers,ATcontent,GCcontent,Basecount,Transcribe\n""RNA:Seqtype,Kmers,ATcontent,GCcontent,Basecount,Translate\n""Protein:Seqtype,Kmers,Hydrophobicity,MolecularWeight)\n")

    #Based on the sequence and method user inputs the relevant operation is done
    #get sequence type
    if Operation == 'Seqtype':
        Type = sequence_analyzer.Sequence.getSeqType(sequence)
        print('\n Type of the sequence is:',Type,'\n')

    #get Kmers in a sequence
    elif Operation == 'Kmers':
        K = int(input('input K value:\n'))
        Kmers = object.getKmers(K)
        print('\n','Kmers of the sequence are:','\n',Kmers,'\n')

    #get AT content of a sequence
    elif Operation == 'ATcontent':
        ATcontent = object.ATcontent
        print('\n','AT content of the sequence:','\n',ATcontent,'\n')

    #get GC content of a sequence
    elif Operation == 'GCcontent':
        GCcontent = object.GCcontent
        print('\n','GC content of the sequence:','\n',GCcontent,'\n')

    #get base count of a sequence
    elif Operation == 'Basecount':
        Basecount = object.Basecount
        print('\n','Count of each base of the sequence:','\n',Basecount,'\n','A graph comparing each base is generated')

    #get Transcribed sequence
    elif Operation == 'Transcribe':
        Transcribe = object.Transcribed_seq
        print('\n','Transcribed sequence is','\n',Transcribe,'\n')
        #the transcribed RNA sequence can be used as the input sequence and operations can be done if the user prefers
        seq = input('Find details about transcribed RNA sequence?(Yes/No)\n')
        if seq == 'Yes':
            RNA = sequence_analyzer.RNAseq(Transcribe)
            Seqtype = sequence_analyzer.Sequence.getSeqType(Transcribe)
            #get details of the trasncribed RNA sequence
            AUcontent = RNA.ATcontent
            GCcontent = RNA.GCcontent
            Basecount = RNA.Basecount
            Translate = RNA.TranslatedSeq
            print('Sequence type is:',Seqtype)
            print('AU content is:',AUcontent)
            print('GC content is:',GCcontent)
            print('Count of each base is:',Basecount)
            print('Translated sequence is:',Translate)
        #if the user does not want to check details of the transcribed sequence program ends
        if seq == 'No':
            print("\nEnding the program")

    #get translated sequence
    elif Operation == 'Translate':
        Translate = object.TranslatedSeq
        print('\n Translated sequence is','\n',Translate,'\n')

    #get Hydrophobicity of a sequence
    elif Operation == 'Hydrophobicity':
        Hydrophobicity = object.Hydrophobicity
        print('\n','Hydrophobicity of the sequence is:','\n',Hydrophobicity,'\n')


    #get Molecular weight of a sequence
    elif Operation == 'MolecularWeight':
        Molecularweight = object.MolecularWeight
        print('\n','Molecular weight is:','\n',Molecularweight,'\n')

    #If the user inputs an unavailable method or inputs the method with incorrect spellings,cases a error message will appear
    else:
        print('Check the spellings and cases of the input method')

#if the user inputs two files comparsions of AT and GC contents can be done
else:
    file_1,file_2 = files
    print('You entered two files:',file_1,file_2)

    sequence1 = sequence_analyzer.Sequence.fastasplit(file_1)
    sequence2 = sequence_analyzer.Sequence.fastasplit(file_2)

    sequence_1_type = sequence_analyzer.Sequence.getSeqType(sequence1)
    sequence_2_type = sequence_analyzer.Sequence.getSeqType(sequence2)

    #only sequences of the same type can be compared
    if sequence_1_type == 'DNA' and sequence_2_type == 'DNA':
        comparison = sequence_analyzer.DNAsequence.compareATandGCcontents(sequence1, sequence2)
        print('\nComparison graph  of AT and GC contents is created')

    elif sequence_1_type == 'RNA' and sequence_2_type == 'RNA':
        comparison = sequence_analyzer.RNAseq.compareATandGCcontents(sequence1, sequence2)
        print('\nComparison graph of AT and GC contents is created')

    #if the sequences are differnt an error message will appear
    else:
        print('This program cannot compare a RNA and DNA sequence')

print('----------------------------------------------------------------------------------------------------------------')