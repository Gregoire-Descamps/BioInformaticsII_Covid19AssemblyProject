from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt


def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    for base in sequence:
        if base not in complement_dict:
            raise ValueError('Error, sequence includes at least an invalid character {}'.format({base}))
    reverse_sequence =sequence[::-1]
    reverse_complement_sequence=''.join(complement_dict.get(base.upper(), base.upper()) for base in reverse_sequence)
    return reverse_complement_sequence


def LoadSeqIOFast(file, format):
    # load a fastq file into an iterable SeqIO object
    obj = SeqIO.parse(file, format)
    return obj


def kmerDict_creation_from_fasta(file, kmer_length):
    SeqIO = LoadSeqIOFast(file, 'fasta')
    kmer_dict = {}
    for seq in SeqIO:
        maxKmer = len(seq.seq) - kmer_length + 1
        if len(seq.seq) < kmer_length:
            raise Exception(NameError('sequence length in file too short'), f'{seq.id} sequence is too short')

        for i in range(maxKmer):
            kmer = str(seq.seq[i:i + kmer_length]).upper()
            rKmer = reverse_complement(kmer)

            if kmer in kmer_dict:
                kmer_dict[kmer] += 1

            elif rKmer in kmer_dict:
                kmer_dict[rKmer] += 1

            else:
                kmer_dict[kmer] = 1

    return kmer_dict


def kmerDict_creation_from_fastq(file, kmer_length):
    SeqIO = LoadSeqIOFast(file, 'fastq')
    kmer_dict = {}

    for seq in SeqIO:
        maxKmer = len(seq.seq) - kmer_length + 1
        if len(seq.seq) < kmer_length:
            raise Exception(NameError('sequence length in file too short'), f'{seq.id} sequence is too short')

        for i in range(maxKmer):
            kmer = str(seq.seq[i:i + kmer_length]).upper()
            rKmer = reverse_complement(kmer)

            if kmer in kmer_dict:
                kmer_dict[kmer] += 1

            elif rKmer in kmer_dict:
                # Only 1 rkmer or kmer are kept, Kmers are considered unique within the genome
                # and without reverse complementary sequence in another location.
                kmer_dict[rKmer] += 1

            else:
                kmer_dict[kmer] = 1

    return kmer_dict

def KmerDictMetrics(dict, range, toTitle):

    print(f'number of Kmer in dict : {len(dict)}')
    print(f'Total number of occurrences in dict detected : {sum(dict.values())}')

    plt.hist(dict.values(), bins=range[1]-range[0], range=range)
    plt.title(f'{toTitle} kmers (k={len(list(dict.keys())[0])}) occurrences')
    plt.show()
    print(f'max count for a kmer : {max(dict.values())}')


def DictTrimming(dict):
    # instead of removing them, under-threshold Kmer are set to False
    clearDict = {}
    for key, value in dict.items():
        if value < 20:
            clearDict[key] = False
        else:
            clearDict[key] = value

    return clearDict


def KmerContentCompare(dictQuery, dictSubject):
    if len(list(dictQuery.keys())[0]) != len(list(dictQuery.keys())[0]):
        raise Exception(f"Length of the Kmer is the 2 dictionaries to compare are different! {len(list(dictQuery.keys())[0])} VS {len(list(dictSubject.keys())[0])} ")
    indict = 0
    notindict = 0

    for key, value in dictSubject.items():
        if key in dictQuery and value != False:
            indict += 1
        elif reverse_complement(key) in dictQuery and value != False:
            indict += 1
        elif value:
            notindict += 1

    total = indict + notindict
    print(f'\nfound {indict} kmer in both dict and {notindict} kmer from dict1 are not found in dict2 ({(indict / total * 100)}% accuracy)')
    print(f'{len(dictQuery)-indict} kmer are missing in second dict ({indict/len(dictQuery)*100}% coverage)\n')


def branchScoring(branch, KmerDict):
    score = 0
    for Kmer in branch:
        if Kmer in KmerDict:
            score += KmerDict[Kmer]
        elif reverse_complement(Kmer) in KmerDict:
            score += KmerDict[reverse_complement(Kmer)]
        else:
            raise Exception(f'Error, Kmer {Kmer} is not in Dict!')

    return score / len(branch)


# NewNewGraphConstructor
def ContigConstructor(inpKmerDict, outputFileName='output.fasta'):
    temp_KmerDict = inpKmerDict.copy()
    contigCount = 0
    contigs = []
    bases = ['A', 'T', 'G', 'C']
    # creating or emptying the output file
    outputfile = open(outputFileName, 'w')
    outputfile.write('')
    outputfile.close()
    outputfile = open(outputFileName, 'a')

    for Kmer, occurences in temp_KmerDict.items():
        # Seed only on k-mer with at least 60 occurrences
        if occurences != False and occurences >= 60:
            temp_KmerDict[Kmer] = False
            branches = []
            firstSeed = []
            firstSeed.append(Kmer)
            branches.append(firstSeed)
            continueElong = True
            completedBranches = []
            KmerSet = set()

            while continueElong:
                continueElong = False
                tempBranches = []
                index = 0
                completedIndexes = []

                for branch in branches:
                    endKmerList = []
                    starKmerList = []
                    lastKmer = branch[len(branch) - 1]
                    firstKmer = branch[0]
                    KmerToAppend = []
                    KmerToInsert = []
                    found = False

                    for base in bases:
                        endKmerList.append(lastKmer[1:len(lastKmer)] + base)
                        endKmerList.append(reverse_complement(endKmerList[-1]))
                        starKmerList.append(base + firstKmer[0:len(firstKmer) - 1])
                        starKmerList.append(reverse_complement(starKmerList[-1]))

                    for i in range(1, 5):
                        # elongating kmer's 5' end
                        if endKmerList[2 * i - 2] in temp_KmerDict and temp_KmerDict[endKmerList[2 * i - 2]] and not (endKmerList[2 * i - 2] in KmerSet):
                            continueElong = found = True
                            KmerToAppend.append(endKmerList[2 * i - 2])
                            KmerSet.add(endKmerList[2 * i - 2])


                        elif endKmerList[2 * i - 1] in temp_KmerDict and temp_KmerDict[endKmerList[2 * i - 1]] and not (endKmerList[2 * i - 1] in KmerSet):
                            continueElong = found = True
                            KmerToAppend.append(endKmerList[2 * i - 2])
                            KmerSet.add(endKmerList[2 * i - 1])

                        # elongating Kmer's 5' end
                        if starKmerList[2 * i - 2] in temp_KmerDict and temp_KmerDict[starKmerList[2 * i - 2]] and not (starKmerList[2 * i - 2] in KmerSet):
                            continueElong = found = True
                            KmerToInsert.append(starKmerList[2 * i - 2])
                            KmerSet.add(starKmerList[2 * i - 2])


                        elif starKmerList[2 * i - 1] in temp_KmerDict and temp_KmerDict[starKmerList[2 * i - 1]] and not (starKmerList[2 * i - 1] in KmerSet):
                            continueElong = found = True
                            KmerToInsert.append(starKmerList[2 * i - 2])
                            KmerSet.add(starKmerList[2 * i - 1])

                    # not finding a Kmer signify not elongated branch (full branch)
                    # full branches indexes are kept to remove them after iteration
                    if not found:
                        completedIndexes.append(index)

                    # scoring found corresponding Kmer and setting them all to false
                    maxScore = 0
                    maxIndex = None
                    index = 0
                    for Kmer in KmerToInsert:

                        if Kmer in temp_KmerDict:
                            if temp_KmerDict[Kmer] > maxScore:
                                maxScore = temp_KmerDict[Kmer]
                                maxIndex = index
                            temp_KmerDict[Kmer] = False

                        elif reverse_complement(Kmer) in temp_KmerDict:
                            if temp_KmerDict[reverse_complement(Kmer)] > maxScore:
                                maxScore = temp_KmerDict[reverse_complement(Kmer)]
                                maxIndex = index
                            temp_KmerDict[reverse_complement(Kmer)] = False

                        index += 1

                    # adding max Kmer to the main branch, other serves as seed for new branches
                    for i in range(len(KmerToInsert)):
                        if i == maxIndex:
                            branch.insert(0, KmerToInsert[i])
                        else:
                            tempBranches.append([KmerToInsert[i]])

                    # scoring found corresponding Kmer and setting them all to false
                    maxScore = 0
                    maxIndex = None
                    index = 0
                    for Kmer in KmerToAppend:

                        if Kmer in temp_KmerDict:
                            if temp_KmerDict[Kmer] > maxScore:
                                maxScore = temp_KmerDict[Kmer]
                                maxIndex = index
                            temp_KmerDict[Kmer] = False

                        elif reverse_complement(Kmer) in temp_KmerDict:
                            if temp_KmerDict[reverse_complement(Kmer)] > maxScore:
                                maxScore = temp_KmerDict[reverse_complement(Kmer)]
                                maxIndex = index
                            temp_KmerDict[reverse_complement(Kmer)] = False

                        index += 1

                    # adding max Kmer to the main branch, other serves as seed for new branches
                    # using a temporary list for new branches as it is impossible/dangerous to modify a list when iterating on it
                    for i in range(len(KmerToAppend)):
                        if i == maxIndex:
                            branch.append(KmerToAppend[i])
                        else:
                            tempBranches.append([KmerToAppend[i]])

                if len(tempBranches) > 0:
                    branches += tempBranches

                # removing complete branches from active list to reduce computing time
                if len(completedIndexes) > 0:
                    completedIndexes.sort(reverse=True)
                    for index in completedIndexes:
                        completedBranches.append(branches[index])
                        branches.pop(index)
                    completedIndexes = []

            # adding back complete branches
            branches += completedBranches

            # Construct contigs
            for branch in branches:
                contigList = branch
                branchScore = branchScoring(contigList, inpKmerDict)


                # ignonring contigs shorter than (2 kmer) length  and a score below 500
                if len(contigList) >= 2 * (len(contigList[0])) + 1 and branchScore > 400:
                    contigCount += 1
                    seqContig = contigList[0]

                    for i in range(1, len(contigList)):
                        seqContig += contigList[i][-1]
                    contigs.append(seqContig)
                    contig = SeqIO.SeqRecord(Seq(seqContig), 'Contig.' + str(contigCount), '', 'Severe acute respiratory syndrome coronavirus 2 isolate SARS-CoV-2/human/VNM/nCoV-19-01S/2020 Debruijn graph contigs')
                    SeqIO.write(contig, outputfile, 'fasta')

                # "Releasing" highscore Kmer from disrupted contigs by re-adding their score
                else:
                    for Kmer in contigList:
                        if Kmer in temp_KmerDict:
                            if inpKmerDict[Kmer] > 100:
                                temp_KmerDict[Kmer] = inpKmerDict[Kmer]
                        else:
                            if inpKmerDict[reverse_complement(Kmer)] > 100:
                                temp_KmerDict[reverse_complement(Kmer)] = inpKmerDict[reverse_complement(Kmer)]


    return contigs


def AssemblyQualityAssess(fastafile, assemblyDict):
    DictCompare = kmerDict_creation_from_fasta(fastafile, len(list(assemblyDict.keys())[0]))
    KmerContentCompare(assemblyDict, DictCompare)

    print('Metrics for the assembly dict')
    KmerDictMetrics(DictCompare, [0, 100],str(fastafile)+' ')


def StatsCalculator(assemblyFile, genomeFile):
    contigs = []
    genomefasta = LoadSeqIOFast(genomeFile, 'fasta')
    assemblyfasta = LoadSeqIOFast(assemblyFile, 'fasta')

    for seq in genomefasta:
        genomeSize = len(seq.seq)

    for seq in assemblyfasta:
        contigs.append(seq.seq)

    contigs.sort(reverse=True, key=len)
    covPercent = 0
    contignum = 1
    while covPercent < 50:
        covPercent += len(contigs[contignum-1]) / genomeSize * 100
        N50 = len(contigs[contignum - 1])
        L50 = contignum
        contignum+=1
    print(f'Assembly N50 : {N50}\nAssembly L50 : {L50}')


# Merge contigs in fasta file if they have common Kmer (overlap > Kmer len)
# recreate a dictionary of K length from the assembly fasta file
# Kmer are scored by using the occurrences they are found time 500 (above thresholds for seeding and releasing Kmer, and above scoring mean)
def Contigmerging(intputFasta, lenKmer, outputFile):
    assemblyKmerDict = kmerDict_creation_from_fasta(intputFasta, lenKmer)

    for key in assemblyKmerDict.keys():
        assemblyKmerDict[key] *= 500

    return ContigConstructor(assemblyKmerDict, outputFile)



readfileName ='SRR21719088_1.fastq'


# Building Dictionaries of Kmer of K = 20 ,25 30
# Then printing frequency of Kmer
print('\nDict20:\n')
print('\t----------')
Dict20 = kmerDict_creation_from_fastq(readfileName, 20)
KmerDictMetrics(Dict20,[1,100], "Normal")
print('\nDict25:\n')
print('\t----------')
Dict25 = kmerDict_creation_from_fastq(readfileName, 25)
KmerDictMetrics(Dict25,[1,100],"Normal")
print('\nDict30:\n')
print('\t----------')
Dict30= kmerDict_creation_from_fastq(readfileName, 30)
KmerDictMetrics(Dict30,[1,100],"Normal")


# Cleaning Dictionaries of Kmer of K = 20 ,25 30 (Threshold = 20)
# Then printing frequency of Kmer
print('\nClearDict20:\n')
print('\t----------')
ClearDict20 = DictTrimming(Dict20)
KmerDictMetrics(ClearDict20, [2,101],"Clear")
print('\nClearDict25:\n')
print('\t----------')
ClearDict25 = DictTrimming(Dict25)
KmerDictMetrics(ClearDict25, [2,101],"Clear")
print('\nClearDict30:\n')
print('\t----------')
ClearDict30 = DictTrimming(Dict30)
KmerDictMetrics(ClearDict30, [2,101],"Clear")

#creating Dictionaries from the SARSCoV2 Genome
print('DictCovid25:\n')
print('\t----------')
DictCovid25 = kmerDict_creation_from_fasta('GCA_011545535.1.fasta', 25)
print('DictCovid30:\n')
print('\t----------')
DictCovid30 = kmerDict_creation_from_fasta('GCA_011545535.1.fasta', 30)


# Compare the different dictionaries with the kmer from the genome
# Accuracy is the percentage of (active) reads Kmer that are in the genome
# Coverage is the percentage of genome's Kmer present (and active) in the read dictionary
print('\nComparing kmer dict 25 with genome kmer dict :\n')
print('\t----------')
KmerContentCompare(DictCovid25, Dict25)
print('\nComparing kmer dict 30 with genome kmer dict :\n')
print('\t----------')
KmerContentCompare(DictCovid30, Dict30)
print('\nComparing Clean kmer dict 25 with genome kmer dict :\n')
print('\t----------')
KmerContentCompare(DictCovid25, ClearDict25)
print('\nComparing Clean kmer dict 30 with genome kmer dict :\n')
print('\t----------')
KmerContentCompare(DictCovid30, ClearDict30)

# Dictionary for quality assessment only, it show the genome's Kmer occurrences in the reads dictionary
CurratedDict30 = {}
revKeyInDict = 0
fwdKeyInDict = 0

for key, value in Dict30.items():
    if key in DictCovid30:
        CurratedDict30[key] = value
        fwdKeyInDict += 1
    elif reverse_complement(key) in DictCovid30:
        CurratedDict30[reverse_complement(key)] = value
        revKeyInDict += 1

print('\n\nMetrics for genome Kmer in the K =30 read dictionary\n')
print('\t----------')
KmerDictMetrics(CurratedDict30,[1,100],"Genome")
KmerDictMetrics(CurratedDict30,[1,2000],"Genome")
KmerContentCompare(DictCovid30,CurratedDict30)
# We can see that beetween 15 and ~75 occurrences, there is no Kmer from genome but about 100 "true" Kmer have only 1 occurrence
# The mode for the occurrences seems to be around 450-500 for the genome's kmer



#assemblies building
print('\nAssembly with k=20')
print('\t----------')
assemblyK20 = ContigConstructor(ClearDict20, 'SARS-CoV-2_Assembly_K20.fasta')
print(f"Number of Contigs in assembly : {len(assemblyK20)}")
AssemblyQualityAssess('SARS-CoV-2_Assembly_K20.fasta',DictCovid30)
StatsCalculator('SARS-CoV-2_Assembly_K20.fasta', 'GCA_011545535.1.fasta')

print('\nAssembly with k=25')
print('\t----------')
assemblyK25 = ContigConstructor(ClearDict25, 'SARS-CoV-2_Assembly_K25.fasta')
print(f"Number of Contigs in assembly : {len(assemblyK25)}")
AssemblyQualityAssess('SARS-CoV-2_Assembly_K25.fasta',DictCovid30)
StatsCalculator('SARS-CoV-2_Assembly_K25.fasta', 'GCA_011545535.1.fasta')

print('\nAssembly with k=30')
print('\t----------')
assemblyK30 = ContigConstructor(ClearDict30, 'SARS-CoV-2_Assembly_K30.fasta')
print(f"Number of Contigs in assembly : {len(assemblyK30)}")
AssemblyQualityAssess('SARS-CoV-2_Assembly_K30.fasta',DictCovid30)
StatsCalculator('SARS-CoV-2_Assembly_K30.fasta', 'GCA_011545535.1.fasta')


# merge the Contigs in the k30 assembly until the number of contigs stop reducing
mergedContigs = Contigmerging('SARS-CoV-2_Assembly_K30.fasta', 25, 'SARS-CoV-2_Assembly_K30_Full.fasta')
lastLength = len(mergedContigs) + 1

while lastLength > len(mergedContigs):
    lastLength = len(mergedContigs)
    mergedContigs = Contigmerging('SARS-CoV-2_Assembly_K30_Full.fasta', 25, 'SARS-CoV-2_Assembly_K30_Full.fasta')
print('Merged contigs Full assembly:')
print('\t----------')
AssemblyQualityAssess('SARS-CoV-2_Assembly_K30_Full.fasta',DictCovid30)
StatsCalculator('SARS-CoV-2_Assembly_K30_Full.fasta', 'GCA_011545535.1.fasta')

