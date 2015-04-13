import gzip
from modtransform.mod.stream import gen_mod

def do_nothing(data, delta):
    return delta

def addition(data, delta):
    return delta + len(data)

def subtraction(data, delta):
    return delta - len(data)

def error_handler(data, delta):
    exit("Not controlling all options")

def build_transform(mod_file, reverse=False):
    """Return nested dictionary: chr: pos1: pos2.
    
    :reverse will map transform in opposite direction, must be type boolean

    """
    assert isinstance(reverse, bool), "reverse type error - boolean required"
    mod = gen_mod(mod_file)
    transform = {}

    if reverse:
        adjustment_direction = {'s': do_nothing, 'd': addition, 'i': subtraction}
    else:
        adjustment_direction = {'s': do_nothing, 'd': subtraction, 'i': addition}
    
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
        except:
            transform[chrom] = [(int(data[2]), delta),]
    return transform

