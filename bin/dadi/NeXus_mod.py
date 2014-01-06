import re
import dadi
reload(dadi.Spectrum_mod)
reload(dadi)

"""
Written by R. Gutenkunst. Modified by M. G. Harvey to accept alignments containing missing data 
(non-ACGT).

Usage: python NeXus_mod.py

"""

def data_dict_from_file(filename, pop_assignments):
    """
    Parse data from NeXus file into data dictionary.

    filename: Name of NeXus file
    pop_assignments: Dictionary mapping taxon identifiers to corresponding
                     populations. If one of the taxa is to be used as the
                     outgroup, it should be assigned to population 'outgroup'.
    """
    # First parse the data from the file
    with open(filename) as fid:
        # Check header to ensure NeXus file
        line = fid.readline()
        if not line.startswith('#NEXUS'):
            raise ValueError('File {0} must be in NeXus format.'.format(filename))

        # Read until the 'Begin data' line
        line = fid.readline()
        while not line.startswith('Begin data'):
            line = fid.readline()

        # On the next Dimensions line, extract the number of taxa and characters
        line = fid.readline()
        ntax = int(re.search(r'ntax=(\d*)', line).group(1))
        nchar = int(re.search(r'nchar=(\d*)', line).group(1))

        # Skip next two lines
        line = fid.readline()
        line = fid.readline()

        # Run through lines with taxon data, copying data into taxa and
        # sequencees lists.
        taxa,sequences = [],[]
        for t in range(ntax):
            line = fid.readline()
            taxon, seq = line.split()
            seq = seq.upper()
            taxa.append(taxon)
            sequences.append(seq)

    # Count number of populations in data
    populations = set(pop_assignments.values())
    populations.discard('outgroup')
    
    # Determine which taxon is the outgroup.
    for taxon,pop in pop_assignments.items():
        if pop=='outgroup':
            outgroup_ii = taxa.index(taxon)
            break
    else:
        outgroup_ii = None
    
    # Iterate through sites in the alignment, creating a data_dict entry
    # for each.
    data_dict = {}
    for site in range(nchar):
        # Pull out all the calls for that site
        calls = [s[site] for s in sequences]

        # Determine what alleles are present at that site 
        palleles = set(calls).intersection('ACTG-N?')
    
        site_name = 'site_{0}'.format(site+1)

        alleles = set(calls).intersection('ACTG')
    
        if len(alleles) > 2:
            print('Ignoring site {0}, which has more than two '
                  'alleles'.format(site+1))
            continue
        elif len(alleles) < 2:
            print('Ignoring site {0}, which has less than two '
                  'alleles'.format(site+1))
            continue

        entry = {}
        data_dict[site_name] = entry

        # Record which alleles are segregating.
        segregating = tuple(alleles)
        entry['segregating'] = segregating
    
        # Record outgroup allele
        if outgroup_ii is not None:
            entry['outgroup_allele'] = calls[outgroup_ii]

        # Initialize dictionary that will hold call counts for each population.
        call_dict = {}
        entry['calls'] = call_dict
        for pop in populations:
            call_dict[pop] = [0,0]
    
        # Run through the calls, assigning each the the proper population and
        # allele.
        for taxon, call in zip(taxa, calls):
            pop = pop_assignments[taxon]
            # Ignore the outgroup here
            if pop is 'outgroup':
                continue
            if call in segregating:
                if call == segregating[0]:
                    aa = 0
                elif call == segregating[1]:
                    aa = 1
                call_dict[pop][aa] += 1

    return data_dict
