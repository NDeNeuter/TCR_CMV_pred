#!/usr/bin/env python3

from collections import Counter

#############
## CLASSES ##
#############

class TCR_Repertoire():
    
    """ Class representing a collection of T-cells.
    Takes as input:
    - TCR_list: a list of TCR instances
    - name: use to identify repertoire """
    
    def __init__(self, TCR_list, name = ''):
        
        self.name = name
        self.TCRs = TCR_list
        self.unique_cdr3_sequences = []
        
    
    def __len__(self):
        
        """ Number of unique TCRs in the repertoire.
        Some TCRs might be present as seperate TCR instances if they're encoded by different DNA sequences. """
        
        return len(set(self.TCRs))
    
    
    def add_TCR(self, tcr_instance):
        
        """ Add a TCR instance to the repertoire. """
        
        self.TCRs.append(tcr_instance)
        
        
    def unique_cdr3_counts(self):
        
        """ Return a dict of unique CDR3 sequences (keys) and how often they occur in the repertoire (value). """
        
        cdr3_dict = {}
        for tcr in self.TCRs:
            cdr3_dict.setdefault(tcr.cdr3, 0)
            cdr3_dict[tcr.cdr3] += 1
    
        return cdr3_dict
    
    
    def unique_cdr3_read_count(self):
        
        """ Return a dict of unique CDR3 sequences (keys) and their total read count as observed in the repertoire (value). """
        
        cdr3_read_dict = {}
        for tcr in self.TCRs:
            cdr3_read_dict.setdefault(tcr.cdr3, 0)
            cdr3_read_dict[tcr.cdr3] += tcr.read_count
            
        return cdr3_read_dict
    
    
    def total_read_count(self):
        
        """ Returns the total number of reads in the repertoire. """
        
        return sum([tcr.read_count for tcr in self.TCRs])
    
    
    def make_productive_repertoire(self):
        
        """ Limit the TCRs in the repertoire to those with productive TCR sequences
        and return them as a new repertoire instance. """
        
        productive_list = []
        for tcr in self.TCRs:
            if tcr.frame == 'In':
                productive_list.append(tcr)
                
        return TCR_Repertoire(productive_list, '{}_productive'.format(self.name))
        
        
    def write_as_transactions(self, filepath):
        
        """ Write all TCRs in repertoire to a file.
        Format TCRs as transaction lists. """
        
        with open(filepath, 'w') as out:
            
            header = '# TCR transaction data'
            out.write(header)
            for tcr in self.TCRs:
                out.write('\n'+tcr.as_transaction())
        
        
class TCR():
    
    """ Class representing T-cell receptors based on data retrieved from the ImmuneAccess database.
    Takes as input:
    - dataline: a single (non-header) line from an ImmuneAccess datafile
    - filepath: path to the file containing the data (required for meta info)
    - format: not all files from ImmuneAccess use the same format, use to specify which kind of format the file is in.
    If using a new format -> add it to code in __init__ method. """
    
    def __init__(self, datadict, file_origin=None):

        self.data = datadict
        self.filename = file_origin
        self.nt_sequence = datadict['nucleotides']
        self.cdr3 = datadict['cdr3']
        self.read_count = int(datadict['read_count'])
        self.frame = datadict['frame']
        self.v_max = datadict['v_max']
        self.v_family = datadict['v_family']
        self.v_gene = datadict['v_gene']
        self.j_max = datadict['j_max']
        self.j_family = datadict['j_family']
        self.j_gene = datadict['j_gene']

        # set CD based on filename (can be either CD4 or CD8)
        self.cd = None
        if 'CD4' in self.filename:
            self.cd = 'CD4'
        elif 'CD8' in self.filename:
            self.cd = 'CD8'

        # set t-cell type based on filename (can be either naive or memory t-cells)
        self.t_cell_type = None
        if 'NAIVE' in self.filename or 'Naive' in self.filename:
            self.t_cell_type = 'naive'
        elif 'TEM' in self.filename or 'Memory' in self.filename:
            self.t_cell_type = 'memory'
        else:
            self.t_cell_type = ''
  
    
    def __hash__(self):
        
        return hash((self.cdr3, self.v_max, self.j_max))


    def __eq__(self, other):

        return '/'.join([self.v_max, self.cdr3, self.j_max]) == other.__repr__()
        #return self.cdr3 == other.cdr3 and self.v_max == other.v_max and self.j_max == other.j_max

    
    def __ne__(self, other):
    
        return not self.__eq__(other)
    
    
    def __str__(self):
        
        return self.__repr__()
        #return '<TCR instance> - {}, {}, {}'.format(self.v_max, self.cdr3, self.j_max)
    
    
    def __repr__(self):
        
        return '/'.join([self.v_max, self.cdr3, self.j_max])
    
    
    def as_transaction(self, k_mer_length = 3):
        
        """ Return a string of TCR properties (=items) that can be used as transaction data for FIM.
        Items include:
        - consecutive amino acid k-mers of length k_mer_length from the TCR CDR3 sequence
        - V/J gene/family.
        - cell type
        - CD molecule (CD4+ or CD8+) """
        
        cdr3_items = [self.cdr3[x:x+k_mer_length] for x in range(len(self.cdr3)-k_mer_length)]
        transaction_items = [self.cd, self.t_cell_type, self.v_family, self.v_gene, self.j_family, self.j_gene] + cdr3_items
        
        return ','.join([x for x in transaction_items if x != None])
    
        
###############
## FUNCTIONS ##
###############


def read_immuneaccess_data(file_path, version='v1'):
    
    """ Read in TCR data downloaded from ImmuneACCESS and return a list of TCR instances.
    Arguments:
    - file_path: path to the file
    - version: version chosen when downloading data from ImmuneACCESS (v1 or v2).
    Returns:
    - a TCR_Repertoire instance containing all the t-cell receptors in the file
    """
                
    tcrs = []
    
    with open(file_path, 'r') as f:
        
        header = f.readline().strip().split('\t')
        data = f.read().strip().split('\n')
        
        if version == 'v1':
            items = {'rearrangement': 'nucleotides',
                     'amino_acid': 'cdr3',
                     'frame_type': 'frame',
                     'templates': 'read_count',
                     'v_family': 'v_family',
                     'v_gene': 'v_gene',
                     'v_allele': 'v_allele',
                     'j_family': 'j_family',
                     'j_gene': 'j_gene',
                     'j_allele': 'j_allele'}
        elif version == 'v2':
            items = {'nucleotide': 'nucleotides',
                     'aminoAcid': 'cdr3',
                     'sequenceStatus': 'frame',
                     'count (templates/reads)': 'read_count',
                     'vMaxResolved': 'v_max',
                     'vFamilyName': 'v_family',
                     'vGeneName': 'v_gene',
                     'jMaxResolved': 'j_max',
                     'jFamilyName': 'j_family',
                     'jGeneName': 'j_gene'}
        
        header_positions = {}
        for item in items.keys():
            header_positions[items[item]] = header.index(item)
        
        for line in data:
            line = line.strip().split('\t')
            tcr_data = {}
            for item in header_positions.keys():
                tcr_data[item] = line[header_positions[item]]
            if version == 'v1':
                tcr_data['v_max'] = make_max_allele(tcr_data, gene='v')
                tcr_data['j_max'] = make_max_allele(tcr_data, gene='j')
            tcrs.append(TCR(tcr_data, file_origin=file_path))
            
    return TCR_Repertoire(tcrs, name=file_path)

             
def filter_cdr3_sequences(repertoire, cdr3_list):
    
    """ Filter out TCRs from repertoire if their CDR3 sequence is in cdr3_list.
    Returns a new repertoire with the remaining TCR. """
    
    seq_id = {}

    for tcr in repertoire.TCRs:
        seq_id.setdefault(tcr.cdr3, []).append(tcr)
        
    for cdr3 in cdr3_list:
        del seq_id[cdr3]
    
    tcr_list = []
    for value in seq_id.values():
        tcr_list += value
    
    return TCR_Repertoire(tcr_list, name = '{}_non_shared'.format(repertoire.name))


def remove_cdr3_overlap(repertoire_A, repertoire_B):
    
    """ Requires two TCR repertoires as input and removes TCRs from both repertoires if their CDR3 sequence is found in both.
    Returns filtered version of repertoire_A, filtered version of repertoire_B and a list of shared cdr3 sequences"""
    
    print('Total amount of TCRs in each repertoire:\n{}: {}\n{}: {}'\
          .format(repertoire_A.name, len(repertoire_A), repertoire_B.name, len(repertoire_B)))
    
    # get cdr3 sequences from each repertoire
    unique_cdr3_A = repertoire_A.get_unique_cdr3_sequences()
    unique_cdr3_B = repertoire_B.get_unique_cdr3_sequences()
    
    print('Number of unique CDR3 sequences in each repertoire:\n{}: {}\n{}: {}'\
          .format(repertoire_A.name, len(unique_cdr3_A), repertoire_B.name, len(unique_cdr3_B)))
    
    overlapping_cdr3_seqs = [cdr3 for cdr3 in unique_cdr3_A if cdr3 in unique_cdr3_B]
    
    print('Number of unique, shared CDR3 sequences between both repertoires: {}'.format(len(set(overlapping_cdr3_seqs))))

    print('Filtering {}'.format(repertoire_A.name))
    filtered_A = filter_cdr3_sequences(repertoire_A, set(overlapping_cdr3_seqs))
    
    print('Filtering {}'.format(repertoire_B.name))
    filtered_B = filter_cdr3_sequences(repertoire_B, set(overlapping_cdr3_seqs))
    
    return filtered_A, filtered_B, set(overlapping_cdr3_seqs)


def make_max_allele(tcr_data, gene='v'):
    
    """ Used when reading TCR data from v1 ImmuneAccess files.
    Takes information from tcr_data and makes the TCR attribute v_max or j_max (depending on gene).
    Args:
    - tcr_data: dictionary with TCR data
    - gene: 'v' or 'j' gene segment
    """
    
    allele = gene+'_allele'
    if tcr_data[allele]:
        return tcr_data[gene+'_gene']+'*'+tcr_data[allele]
    elif tcr_data[gene+'_gene']:
        return tcr_data[gene+'_gene']
    else:
        return tcr_data[gene+'_family']
