from collections import Counter
from Bio.Data import CodonTable
from typing import TextIO, List
# Create a class for sequence record so that we can quickly assesse sequence ID and sequence
class Record:
    def __init__(self, record):
        self.id = record['name']
        self.seq = record['seq']

def parse_fasta(f: TextIO) -> List[Record]:
    # Create an empty list for storing records
    records = []
    # Set an default empty record
    record = {
        'name': '',
        'seq': ''
    }

    # Set default first line is not the end of record
    end_of_record = False
    for line in f.readlines():
        # Read line by line as while loop runs, and also strip off /n at the end of each line
        line = line.rstrip()
        if not line: # If line is empty, that is not False, skip this loop
            continue

        if line.startswith('>'):
            # If we find the '>' again, that is the end of previous record, so we append the record, and then empty the record
            if end_of_record:
                r = Record(record)
                records.append(r)
                record = {
                    'name': '',
                    'seq': ''
                }
            # Set the end_of_record as True after we pass the first record name
            end_of_record = True
            record['name'] = line.strip('>')

        else:
            # line not startswith '>' are all sequences
            record['seq'] += line
    # Append the last record after the end of while loop
    r = Record(record)
    records.append(r)
    return records

def gc(record: Record):
    # Retrun GC content of record seqeunce, the number is expressed as % and rounded up to 6th decimal number
    nb_counts = Counter(record.seq.upper())
    return round((nb_counts['C'] + nb_counts['G']) / len(record.seq) * 100, 6)

def max_gc(records: List[Record]):
    # Guess the first record has the largest GC content
    max_record = records[0]
    # Compare GC contents of records one by one
    for record in records:
        if gc(record) > gc(max_record):
            max_record = record
    print(max_record.id)
    print(gc(max_record))
    return max_record

def translation(s: str) -> str:
    # Use codon table from biopython
    table = CodonTable.standard_dna_table.forward_table
    stop_codons = CodonTable.standard_dna_table.stop_codons
    protein_seq = ''
    for i in range(0, len(s), 3):
        if len(s[i:i+3]) < 3:
            return 
        if s[i:i+3] in stop_codons:
            return protein_seq
        protein_seq += table[s[i:i+3]]
        i += 3
    

def rc_DNA(s: str) -> str:
    comp_table = {
    'A': 'T',
    'C': 'G',
    'T': 'A',
    'G': 'C'} # Create a complementary table
    rc_s = ''.join([comp_table[i] for i in s])[::-1] # Reverse and complement the DNA sequence
    return rc_s
