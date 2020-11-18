import re
import sys

# ------------------------------ FUNCTIONS --------------------------------- #

def get_fiveutr_coor(gtf_file, chrom):
    """
    Look for the starting and ending position of the 5' UTR in the transcripts
    from protein coding genes belonging to chromosome chrom, given a GTF file.
    Returns a list of these coordinates as tuples, one for each UTR
    """

    # Initialize empty variables
    coor = []
    utrs = []
    prev_gene = ""

    # Read file line by line
    with open(gtf_file, 'r') as file:
        for line in file:

            # The chromosome field can be specified as a number x or as chrx.
            # Unify both cases
            chr = line.split()[0]
            if chr[0:3] == 'chr':
                chr = chr[3:]

            #print(line.split())

            # Extract info from the desired chromosome. We are only interested
            # in the UTR regions from protein coding genes.
            if chr == chrom and line.split()[2] == 'UTR' \
                and line.split()[13] == '"protein_coding";':

                    # Extract relevant fields from the text line
                    this_gene = str(line.split()[9])
                    utr_pos = (int(line.split()[3]), int(line.split()[4]))
                    strand = line.split()[6]

                    # If this UTR belong to the same previous gene, add to list
                    if this_gene == prev_gene:
                        utrs.append(utr_pos)

                    # Else, look for the 5' UTR and save the coordinates
                    else:
                        if utrs:

                            # If UTR coming from + strand, 5'UTR will be the first
                            # in the sequence
                            if strand == '+':
                                five_utr = min(utrs, key=lambda x:x[0])

                            # If UTR coming from - strand, 5'UTR will be the last
                            # in the sequence
                            elif strand == '-':
                                five_utr = max(utrs, key=lambda x:x[0])

                            coor.append(five_utr)

                        # Change to next gene
                        prev_gene = this_gene
                        utrs = []

    return coor



def parse_dna(genome_file, query):
    """
    Given a fasta file, extracts the sequence information from a specific
    given query. Returns a string with the sequence.
    """

    # Don't extract the text in the beginning
    extract = False
    seq = ""

    # Read file line by line
    with open(genome_file, 'r') as file:
        for line in file:

            # Extract line with the sequence. Exclude trailing characters
            if extract:
                seq += str(line.rstrip())

            # Stop extraction after new sequence entry (new chromosome)
            if line[0] == '>':
                extract = False

            # Start extraction after the query line.
            if line[:6] == query:
                extract = True

    return seq


def reverse_complementary(seq):
    """Returns the reverse complementary of a DNA sequence"""

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[n] for n in reversed(seq))



def scan_seq(seq, pattern):
    """
    Given a sequence and a regular expression pattern, search for all
    overlapping occurrences of this pattern in the DNA sequence and its
    reverse complementary sequence
    """

    # Look for matches in the sequence
    matches = [str(match.group(1)) for match in re.finditer(pattern, seq)]

    # Look for matches in the reverse complementary of the sequence
    revcomp_seq = reverse_complementary(seq)
    matches += [str(match.group(1)) for match in re.finditer(pattern, revcomp_seq)]

    return matches



def design_grna(seq):
    """Designs a guide RNA given a DNA template"""

    transcript = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A'}
    grna = "".join(transcript[n] for n in seq)

    return grna



def calculate_hbonds(grna):
    """
    Calculate a score for a guide RNA depending on the number of
    hydrogen bonds that it can form during RNA-DNA hybridization
    """

    scores = {'A': 2, 'U': 2, 'C': 3, 'G': 3}
    v = [scores[n] for n in grna]

    return sum(v)



def predict_targets(dna, grna):
    """
    Predicts and returns all the possible targets of a guide RNA in
    the given DNA sequence.
    """

    # Find targets by the reverse transcription of RNA into cDNA
    revtranscr = {'A': 'T', 'U': 'A', 'G': 'C', 'C': 'G'}
    cdna = "".join(revtranscr[n] for n in grna)
    targets = re.findall(cdna, dna)

    # Find targets in the other DNA strand
    revcomp_cdna = reverse_complementary(cdna)
    targets += re.findall(revcomp_cdna, dna)

    return targets


# ----------------------------------- MAIN -------------------------------------- #
# ---------------------------- Assignment part a) ------------------------------- #

def main():

    # Path to the necessary files
    genome_file = sys.argv[1]
    gtf_file = sys.argv[2]
    chromosome = str(sys.argv[3])

    # Extract DNA sequence from the desired chromosome from the given genome file
    query = '>chr' + chromosome
    dna = parse_dna(genome_file, query)
    print('DNA sequence parsed')
    print(f'The length of the query {query} is {len(dna)} bp')

    # Extract DNA sequence of 5' UTRs of protein coding genes
    utr_coor = get_fiveutr_coor(gtf_file, chromosome)
    print(f'Number of 5\' UTRs found: {len(utr_coor)}')
    utrs = [dna[i - 1:j] for (i, j) in utr_coor]


# ---------------------------- Assignment part b) ------------------------------- #

    # Regex pattern to look for in the UTR sequences for the guide RNA
    pattern = r'(?=([AGCT]{21}GG))'

    # For each 5' UTR (omitting sequences shorter than 22 bp)
    for utr in utrs:
        if len(utr) > 22:

            # Find matches in the 5' UTR to the given pattern
            matches = scan_seq(utr, pattern)

            # Design guide RNAs for the collected 5' UTR DNA matches
            grnas = [design_grna(match[:20]) for match in matches]


# ---------------------------- Assignment part c) ------------------------------- #

            scores = []

            # gRNA ranking for this 5' UTR
            for guide_rna in grnas:

                # Calculate RNA-DNA stability by number of forming hydrogen bonds
                n_hbonds = calculate_hbonds(guide_rna)

                # Predict number of total matches in the chromosome 22 DNA sequence
                total_matches = predict_targets(dna, guide_rna)

                scores.append((n_hbonds - len(total_matches), guide_rna))

            # Print the collected gRNA and its score to the user
            print(sorted(scores, reverse=True))



if __name__ == '__main__':
    main()