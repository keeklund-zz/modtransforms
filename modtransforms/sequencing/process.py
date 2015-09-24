from types import ModuleType
from pysam import AlignmentFile
from modtransforms.formats.bamsam import update_header
from modtransforms.mod.transform import build_transform, find_delta

# will likely move to new location
def get_positions_and_deltas(chrom_mods, curr_chrom, logger):
    """Return tuple containing positional and delta information.

    The first tuple contains positional information.  Every location at which 
    there is a change between species.  The delta tuple is a cumulative measure
    of how much change has occurred up to respective (same index) position.

    """
    assert isinstance(chrom_mods, dict), "chrom_mods must be type dict"
    assert isinstance(curr_chrom, str), "curr_chrom must be type str"
    assert isinstance(logger, ModuleType), "logger must be ModuleType"
    try:
        positions, deltas = zip(*chrom_mods.get(curr_chrom))
        logger.info("CONTIG: '%s'" % curr_chrom)
    except TypeError:
        logger.warn(
            "CONTIG: '%s' is not in MOD File." % \
            curr_chrom)
        positions, deltas = (), ()
    return positions, deltas

def calculate_end_position():
    pass

def atac(args, logger):
    """

    """
    if not args.chrom_sizes:
        exit("Chrom sizes required for bam conversion")

    chrom_mods = build_transform(args.mod, logger)
    input_ = AlignmentFile(args.input, 'rb')

    header = update_header(input_.header.copy(), args.chrom_sizes)
    output = AlignmentFile(args.output, 'wb', header=header)

    curr_chrom = ""
    for line in input_:

        if input_.references[line.reference_id] != curr_chrom:
            curr_chrom = input_.references[line.reference_id]
            positions, deltas = get_positions_and_deltas(chrom_mods,
                                                         curr_chrom,
                                                         logger)
        # if line.is_reverse and (line.reference_length != len(line.seq)):
        #     print line
        #     print line.reference_length
        #     print line.cigar
        #     print len(line.seq)
        #     print len(line.get_reference_positions())
#        try:
        if not line.is_reverse:
            start_delta = find_delta(positions,
                                     deltas,
                                     int(line.reference_start))
            line.reference_start = int(line.reference_start) + start_delta
        else:
            end_delta = find_delta(positions,
                                   deltas,
                                   int(line.reference_end))
            mapped_end = int(line.reference_end) + end_delta
            line.reference_start = mapped_end - len(line.seq) # line.reference_length 
        output.write(line)
#        except IndexError:
#            pass
            
