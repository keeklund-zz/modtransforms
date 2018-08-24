from bisect import bisect_left
from os.path import exists
from types import ModuleType
from modtransforms.mod.stream import gen_mod

def do_nothing(data, delta):
    assert isinstance(data, str), "Data must be type str."
    assert isinstance(delta, int), "Delta must be type int."
    return delta

def addition(data, delta):
    assert isinstance(data, str), "Data must be type str."
    assert isinstance(delta, int), "Delta must be type int."
    return delta + len(data)

def subtraction(data, delta):
    assert isinstance(data, str), "Data must be type str."
    assert isinstance(delta, int), "Delta must be type int."
    return delta - len(data)

def error_handler(data, delta):
    assert isinstance(data, str), "Data must be type str."
    assert isinstance(delta, int), "Delta must be type int."
    exit("Not controlling all options")

def build_transform(mod_file, logger):
    """Return nested dictionary: chr: pos1: pos2.
    
    :mod_file must be type string, and file exist

    :logger is module class 

    """
    assert exists(mod_file), "mod file doesn't exist"
    assert isinstance(logger, ModuleType), "incorrect logger type"
    mod = gen_mod(mod_file)
    transform = {}

    """ No longer used - Karl's original code """
    adjustment_direction = {'s': do_nothing, 'd': addition, 'i': subtraction}
    
    chrom = ''
    delta = 0
    pos = 0
    while mod:
        try:
            data = mod.next().split()
        except:
            break

        if data[1] != chrom:
            delta = 0 
            chrom = data[1]
            """ Initialize beginning of chromosome """
            try:
                transform[chrom].append((1, delta))
            except KeyError:
                transform[chrom] = [(1, delta),]

        """ Deletion in target, so must insert bases to compensate """
        dcount = 0
        curr_pos = int(data[2]) - 1
        while (data[0] == 'd'):
            if dcount == 0:
                pos = int(data[2]) - delta
        
            """ If new chromosome, write out and start next delta """
            if data[1] != chrom: 
                try:
                    transform[chrom].append((pos, delta))
                except KeyError:
                    transform[chrom] = [(pos, delta),]
                delta = 0
                chrom = data[1]
                dcount = 0
                pos = int(data[2]) - delta
                curr_pos = int(data[2]) - 1
                try:
                    transform[chrom].append((1, delta))
                except KeyError:
                    transform[chrom] = [(1, delta),]

            """ Make sure consecutive deletes are for consecutive bases """
            if int(data[2]) != (curr_pos + 1):
                try:
                    transform[chrom].append((pos, delta))
                except KeyError:
                    transform[chrom] = [(pos, delta),]
                dcount = 0
                pos = int(data[2]) - delta

            """ Increment counters for deletion """
            dcount = dcount + 1
            delta = delta + 1
            curr_pos = int(data[2])
            try:
                data = mod.next().split()
            except:
                break                                   

        """ Record deleted bases """
        if (dcount > 0):
            try:
                transform[chrom].append((pos, delta))
            except KeyError:
                transform[chrom] = [(pos, delta),]

        """ Original code from Karl """
        """ handler = adjustment_direction.get(data[0], error_handler) """

        """ Insertion in target, so must delete these inserted bases """
        if (data[0] == 'i'):
            pos = int(data[2]) - delta
            """ Multiple bases can be inserted, and each must be accounted for """
            for i in range(0, len(data[3])): 
                delta = delta - 1
                pos = pos + 1
                try:
                    transform[chrom].append((pos, delta))
                except KeyError:
                    transform[chrom] = [(pos, delta),]

    logger.info("Chromosome MODification transform built")
    return transform

def find_delta(positions, deltas, position):
    """Return accumulated change in position so far on current chromosome.

    """
    assert isinstance(positions, tuple), "positions must be type tuple"
    assert isinstance(deltas, tuple), "deltas must be type tuple"
    assert isinstance(position, int), "position must be type int"
    idx = bisect_left(positions, position)
    try:
        if positions[idx] <= position:
            delta = deltas[idx]
        else:
            delta = deltas[idx - 1]
    except IndexError:
        if len(deltas) > 1:
            delta = deltas[-1]
        else:
            delta = 0
    if idx == 0:
        delta = 0
    return delta
                        
