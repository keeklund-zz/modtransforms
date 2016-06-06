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

    adjustment_direction = {'s': do_nothing, 'd': addition, 'i': subtraction}
    
    chrom = ''
    while mod:
        try:
            data = mod.next().split()
        except:
            break
        if data[1] != chrom:
            delta = 0 
            chrom = data[1]

        handler = adjustment_direction.get(data[0], error_handler)
        delta = handler(data[3], delta)

        try:
            transform[chrom].append((int(data[2]), delta))
        except KeyError:
            transform[chrom] = [(int(data[2]), delta),]
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
                        
