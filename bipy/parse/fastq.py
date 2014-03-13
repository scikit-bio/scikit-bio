__author__ = "Gavin Huttley, Anuj Pahwa"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Anuj Pahwa"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

def MinimalFastqParser(data, strict=True):
    """yields name, seq, qual from fastq file

    Arguments:
        - strict: checks the quality and sequence labels are the same
    """
    if type(data) == str:
        data = open(data)

    # fastq format is very simple, defined by blocks of 4 lines
    line_num = -1
    record = []
    for line in data:
        line_num += 1
        if line_num == 4:
            if strict: # make sure the seq and qual labels match
                assert record[0][1:] == record[2][1:], \
                  'Invalid format: %s -- %s' % (record[0][1:], record[2][1:])
            yield record[0][1:], record[1], record[3]
            
            line_num = 0
            record = []
        
        record.append(line.strip())
    
    if record:
        if strict and record[0]: # make sure the seq and qual labels match
            assert record[0][1:] == record[2][1:], 'Invalid format'
        
        if record[0]: # could be just an empty line at eof
            yield record[0][1:], record[1], record[3]
        
    
    if type(data) == file:
        data.close()

