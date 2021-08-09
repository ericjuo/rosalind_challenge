from collections import Counter
# Create a class for sequence record so that we can quickly assesse sequence ID and sequence
class Record:
    def __init__(self, record):
        self.id = record['name']
        self.seq = record['seq']

def parse_fasta(f):
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

def gc(record):
    # Retrun GC content of record seqeunce, the number is expressed as % and rounded up to 6th decimal number
    nb_counts = Counter(record.seq.upper())
    return round((nb_counts['C'] + nb_counts['G']) / len(record.seq) * 100, 6)

def max_gc(records):
    # Guess the first record has the largest GC content
    max_record = records[0]
    # Compare GC contents of records one by one
    for record in records:
        if gc(record) > gc(max_record):
            max_record = record
    print(max_record.id)
    print(gc(max_record))
    return max_record