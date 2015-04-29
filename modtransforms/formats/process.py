from re import match
from pysam import AlignmentFile
from modtransforms.mod.transform import build_transform, find_delta

def bam(args, logger):
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    curr_chrom = ""

    input_ = AlignmentFile(args.input, 'rb')
    output = AlignmentFile(args.output, 'wb', header=input_.header)

    for line in input_:
        if input_.references[line.reference_id] != curr_chrom:
            curr_chrom = input_.references[line.reference_id]
            try:
                positions, deltas = zip(*chrom_mods.get(curr_chrom))
                logger.info("CONTIG: '%s'" % curr_chrom)
            except TypeError:
                logger.warn(
                    "CONTIG: '%s' is not in MOD File. Skipping." % \
                    curr_chrom)
                positions, deltas = [], 0
                continue
        try:
            start_delta = find_delta(positions,
                                     deltas,
                                     int(line.reference_start))
            end_delta = find_delta(positions,
                                   deltas,
                                   int(line.reference_end))
            if start_delta != end_delta:
                hwm = start_delta
                mod_index = []
                for i,p in enumerate(line.get_reference_positions()):
                    p_delta = find_delta(positions, deltas, p)
                    if hwm != p_delta:
                        mod_index.append((i, hwm-p_delta))
                        hwm = p_delta
                new_cigar = []
                old_cigar = line.cigar

                for mi in mod_index:
                    cigarindex = 0
                    for ci,c in enumerate(old_cigar):
                        
                        if c[1] + mi[1] > 0 and (mi[0] <= cigarindex or mi[0] <= c[1]):
                            new_cigar.append((c[0], c[1] + mi[1]))
                        else:
                            cigarindex = cigarindex + c[1]
                            new_cigar.append(c)
                            
                if len(line.cigar) < len(new_cigar):
                    line.cigar = new_cigar[-1*len(line.cigar):]
                else:
                    line.cigar = new_cigar
            line.reference_start = int(line.reference_start) + start_delta
            output.write(line)
        except IndexError:
            pass

def gtf(args, logger):
    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = 'contig feature source start end score strand frame attributes'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            gene = line.rstrip()
            data = {k:g for k,g in zip(keys.split(), gene.split('\t'))}
            if data.get('contig') != curr_chrom:
                curr_chrom = data.get('contig')
                try:
                    positions, deltas = zip(*chrom_mods.get(curr_chrom))
                    logger.info("CONTIG: '%s'" % curr_chrom)
                except TypeError:
                    logger.warn(
                        "CONTIG: '%s' is not in MOD File. Skipping." % \
                        curr_chrom)
                    positions, deltas = [], 0
                    continue
            try:
                start_delta = find_delta(positions, 
                                         deltas, 
                                         int(data.get('start')))
                end_delta = find_delta(positions, 
                                       deltas,
                                       int(data.get('end')))
                data['start'] = int(data.get('start')) + start_delta
                data['end'] = int(data.get('end')) + end_delta
                out.write('%s\n' % \
                          '\t'.join([str(data.get(k)) for k in keys.split()]))
            except IndexError:
                pass
            
def smrna_12_bed(args, logger):

    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = 'chrom chromStart1 chromEnd1 name score strand chromStart2 chromEnd2 0,0,0 1 0 chromStart3'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            if not line.startswith('chrom'):
                gene = line.rstrip()
                data = {k:g for k,g in zip(keys.split(), gene.split('\t'))}
                if data.get('chrom') != curr_chrom:
                    curr_chrom = data.get('chrom')
                    try:
                        positions, deltas = zip(*chrom_mods.get(curr_chrom))
                        logger.info("CONTIG: '%s'" % curr_chrom)
                    except TypeError:
                        logger.warn(
                            "CONTIG: '%s' is not in MOD File. Skipping." % \
                            curr_chrom)
                        positions, deltas = [], 0
                        continue
                try:
                    start_delta = find_delta(positions, 
                                             deltas, 
                                             int(data.get('chromStart1')))
                    end_delta = find_delta(positions, 
                                           deltas,
                                           int(data.get('chromEnd1')))
                    start_delta = find_delta(positions, 
                                             deltas, 
                                             int(data.get('chromStart2')))
                    end_delta = find_delta(positions, 
                                           deltas,
                                           int(data.get('chromEnd2')))
                    start_delta = find_delta(positions, 
                                             deltas, 
                                             int(data.get('chromStart3')))
                    data['chromStart1'] = int(data.get('chromStart1')) + start_delta
                    data['chromEnd1'] = int(data.get('chromEnd1')) + end_delta
                    data['chromStart2'] = int(data.get('chromStart2')) + start_delta
                    data['chromEnd2'] = int(data.get('chromEnd2')) + end_delta
                    try:
                        data['chromStart3'] = int(data.get('chromStart3')) + start_delta
                    except:
                        pass
                    out.write('%s\n' % \
                              '\t'.join([str(data.get(k)) for k in keys.split()]))
                except IndexError:
                    pass

def smrna_bed(args, logger):

    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = 'chrom chromStart chromEnd name score strand'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            if not line.startswith('chrom'):
                gene = line.rstrip()
                data = {k:g for k,g in zip(keys.split(), gene.split('\t'))}
                if data.get('chrom') != curr_chrom:
                    curr_chrom = data.get('chrom')
                    try:
                        positions, deltas = zip(*chrom_mods.get(curr_chrom))
                        logger.info("CONTIG: '%s'" % curr_chrom)
                    except TypeError:
                        logger.warn(
                            "CONTIG: '%s' is not in MOD File. Skipping." % \
                            curr_chrom)
                        positions, deltas = [], 0
                        continue
                try:
                    start_delta = find_delta(positions, 
                                             deltas, 
                                             int(data.get('chromStart')))
                    end_delta = find_delta(positions, 
                                           deltas,
                                           int(data.get('chromEnd')))
                    data['chromStart'] = int(data.get('chromStart')) + start_delta
                    data['chromEnd'] = int(data.get('chromEnd')) + end_delta
                    out.write('%s\n' % \
                              '\t'.join([str(data.get(k)) for k in keys.split()]))
                except IndexError:
                    pass
            
def smrna_gff3(args, logger):
    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = 'chrom source type start end score strand phase attributes'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            if not line.startswith('#'):
                gene = line.rstrip()
                data = {k:g for k,g in zip(keys.split(), gene.split('\t'))}
                if data.get('chrom') != curr_chrom:
                    curr_chrom = data.get('chrom')
                    try:
                        positions, deltas = zip(*chrom_mods.get(curr_chrom))
                        logger.info("CONTIG: '%s'" % curr_chrom)
                    except TypeError:
                        logger.warn(
                            "CONTIG: '%s' is not in MOD File. Skipping." % \
                            curr_chrom)
                        positions, deltas = [], 0
                        continue
                try:
                    start_delta = find_delta(positions, 
                                             deltas, 
                                             int(data.get('start')))
                    end_delta = find_delta(positions, 
                                           deltas,
                                           int(data.get('end')))
                    data['start'] = int(data.get('start')) + start_delta
                    data['end'] = int(data.get('end')) + end_delta
                    out.write('%s\n' % \
                              '\t'.join([str(data.get(k)) for k in keys.split()]))
                except IndexError:
                    pass
                
def smrna_lib_fa(args, logger):

    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = 'chrom start end strand'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            if line.startswith('>'):
                gene = line.rstrip()

                gene_list = match('>(.*):(\d+)-(\d+)\((.*)\)', gene)
                data = {k:g for k,g in zip(keys.split(), gene_list.groups())}
                if data.get('chrom') != curr_chrom:
                    curr_chrom = data.get('chrom')
                    try:
                        positions, deltas = zip(*chrom_mods.get(curr_chrom))
                        logger.info("CONTIG: '%s'" % curr_chrom)
                    except TypeError:
                        logger.warn(
                            "CONTIG: '%s' is not in MOD File. Skipping." % \
                            curr_chrom)
                        positions, deltas = [], 0
                        continue
                try:
                    start_delta = find_delta(positions, 
                                             deltas, 
                                             int(data.get('start')))
                    end_delta = find_delta(positions, 
                                           deltas,
                                           int(data.get('end')))
                    data['start'] = int(data.get('start')) + start_delta
                    data['end'] = int(data.get('end')) + end_delta
                    out.write('>%s:%s-%s(%s)\n' % \
                              tuple([str(data.get(k)) for k in keys.split()]))
                except IndexError:
                    pass
            else:
                out.write('%s' % line)
                
def smrna_table_txt(args, logger):
    
    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = 'name chrom start end strand mature hairpin'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            gene = line.rstrip()
            data = {k:g for k,g in zip(keys.split(), gene.split())}
            if data.get('chrom') != curr_chrom:
                curr_chrom = data.get('chrom')
                try:
                    positions, deltas = zip(*chrom_mods.get(curr_chrom))
                    logger.info("CONTIG: '%s'" % curr_chrom)
                except TypeError:
                    logger.warn(
                        "CONTIG: '%s' is not in MOD File. Skipping." % \
                        curr_chrom)
                    positions, deltas = [], 0
                    continue
            try:
                start_delta = find_delta(positions, 
                                         deltas, 
                                         int(data.get('start')))
                end_delta = find_delta(positions, 
                                       deltas,
                                       int(data.get('end')))
                data['start'] = int(data.get('start')) + start_delta
                data['end'] = int(data.get('end')) + end_delta
                out.write('%s\n' % \
                          '\t'.join([str(data.get(k)) for k in keys.split()]))
            except IndexError:
                pass

def smrna_txt(args, logger):

    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = 'name num_trna trna_begin trna_end isotype anticodon up_region down_region intron_begin intron_end cove hmm 2_str scanner'
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            gene = line.rstrip()
            data = {k:g for k,g in zip(keys.split(), gene.split())}
            if data.get('name') != curr_chrom:
                curr_chrom = data.get('name')
                try:
                    positions, deltas = zip(*chrom_mods.get(curr_chrom))
                    logger.info("CONTIG: '%s'" % curr_chrom)
                except TypeError:
                    logger.warn(
                        "CONTIG: '%s' is not in MOD File. Skipping." % \
                        curr_chrom)
                    positions, deltas = [], 0
                    continue
            try:
                start_delta = find_delta(positions, 
                                         deltas, 
                                         int(data.get('trna_begin')))
                end_delta = find_delta(positions, 
                                       deltas,
                                       int(data.get('trna_end')))
                data['trna_begin'] = int(data.get('trna_begin')) + start_delta
                data['trna_end'] = int(data.get('trna_end')) + end_delta
                out.write('%s\n' % \
                          '\t'.join([str(data.get(k)) for k in keys.split()]))
            except IndexError:
                pass

