from re import match
from pysam import AlignmentFile
from modtransforms.mod.transform import build_transform, find_delta

def update_header(header, chrom_sizes_file):

    return header

def bam(args, logger):
    if not args.chrom_sizes:
        exit("Chrom sizes required for bam conversion")
        
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    curr_chrom = ""

    input_ = AlignmentFile(args.input, 'rb')

    header = input_.header.copy()
    with open(args.chrom_sizes, "r") as cs_file:
        chrom_sizes = dict(line.split() for line in cs_file)
    for chrom in header['SQ']:
        chrom['LN'] = int(chrom_sizes[chrom['SN']])
    
    output = AlignmentFile(args.output, 'wb', header=header)

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
            # if start_delta != end_delta:
            #     hwm = start_delta
            #     mod_index = []
            #     for i,p in enumerate(line.get_reference_positions()):
            #         p_delta = find_delta(positions, deltas, p)
            #         if hwm != p_delta:
            #             mod_index.append((i, hwm-p_delta))
            #             hwm = p_delta
            #     new_cigar = []
            #     old_cigar = line.cigar

            #     for mi in mod_index:
            #         cigarindex = 0
            #         for ci,c in enumerate(old_cigar):
                        
            #             if len(old_cigar) > 1:
            #                 if c[1] + mi[1] > 0 and (mi[0] <= cigarindex or mi[0] <= c[1]):
            #                     new_cigar.append((c[0], c[1] + mi[1]))
            #                 else:
            #                     cigarindex = cigarindex + c[1]
            #                     new_cigar.append(c)
            #             else:
            #                 new_cigar = old_cigar

            #     if len(line.cigar) < len(new_cigar):
            #         line.cigar = new_cigar[-1*len(line.cigar):]
            #     else:
            #         line.cigar = new_cigar
            line.reference_start = int(line.reference_start) + start_delta
            output.write(line)
        except IndexError:
            pass

def gtf(args, logger):
    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = "seqname source feature start end score strand frame attribute"
    curr_chrom = ""
    with open(args.input, 'r') as input_:
        for line in input_:
            gene = line.rstrip()
            data = gene.split('\t')
            if data[0] != curr_chrom:
                curr_chrom = data[0]
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
                                         int(data[3]))
                end_delta = find_delta(positions, 
                                       deltas,
                                       int(data[4]))
                data[3] = int(data[3]) + start_delta
                data[4] = int(data[4]) + end_delta
                out.write('%s\n' % \
                          '\t'.join([str(d) for d in data]))
            except IndexError:
                pass
            
def smrna_12_bed(args, logger):

    out = open(args.output, 'w')
    chrom_mods = build_transform(args.mod, logger, args.reverse)
    keys = 'chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStart'
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
                    c_start_delta = find_delta(positions, 
                                               deltas, 
                                               int(data.get('chromStart')))
                    c_end_delta = find_delta(positions, 
                                             deltas,
                                             int(data.get('chromEnd')))
                    t_start_delta = find_delta(positions, 
                                               deltas, 
                                               int(data.get('thickStart')))
                    t_end_delta = find_delta(positions, 
                                             deltas,
                                             int(data.get('thickEnd')))
                    b_tmp_data = find_delta(positions, 
                                             deltas, 
                                             (int(data.get('chromStart')) + \
                                              int(data.get('blockStart'))))
                    b_start_delta = b_tmp_data - int(data.get('chromStart'))
                    data['chromStart'] = int(data.get('chromStart')) + c_start_delta
                    data['chromEnd'] = int(data.get('chromEnd')) + c_end_delta
                    data['thickStart'] = int(data.get('thickStart')) + t_start_delta
                    data['thickEnd'] = int(data.get('thickEnd')) + t_end_delta
                    data['blockStart'] = int(data.get('blockStart')) + b_start_delta
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

