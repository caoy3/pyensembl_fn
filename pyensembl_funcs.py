#!/share/public4/data/caoy/softs/python_3.8.12/Python-3.8.12/bin/python3.8

def exon2trans(data,chro,pos,strand,exon_id):
    '''return the list of transcript_ids which contain the exon_id'''
    out = []
    trans_ids = data.transcript_ids_at_locus(chro,pos,strand=strand)
    for i in trans_ids:
        if exon_id in data.exon_ids_of_transcript_id(i):
            out.append(i)
    return out

def get_longest_exon(data,exon_ids):
    '''return the longest exon_id of input exon_ids list'''
    exons_len = []
    for exon_id in exon_ids:
        tmp_exon = data.exon_by_id(exon_id)
        exons_len.append(tmp_exon.length)
    max_idx = [idx for idx,value in enumerate(exons_len) if value == max(exons_len)]
    longest_exons = [exon_ids[idx] for idx in max_idx]
    return longest_exons

def get_trans_length(data,trans_id):
    '''return the total length of exons in the input transcript id'''
    trans = data.transcript_by_id(trans_id)
    exons = trans.exon_intervals
    trans_length = 0
    for exon in exons:
        trans_length += exon[1] - exon[0] + 1
    return trans_length

def get_longest_trans(data,trans_ids):
    '''return the longest trans_id of input trans_ids list'''
    trans_len = []
    for trans_id in trans_ids:
        tmp_trans = data.transcript_by_id(trans_id)
        tmp_exons = tmp_trans.exon_intervals
        tmp_trans_len = 0
        for exon_interval in tmp_exons:
            tmp_trans_len += exon_interval[1] - exon_interval[0] +1
        trans_len.append(tmp_trans_len)
    max_idx = [idx for idx,value in enumerate(trans_len) if value == max(trans_len)]
    longest_trans = [trans_ids[idx] for idx in max_idx]
    return longest_trans

def get_longest_gene(data,gene_ids):
    '''return the longest gene_id of input gene_ids'''
    gene_objs = [data.gene_by_id(gene_id) for gene_id in gene_ids]
    gene_objs_len = [gene_obj.length for gene_obj in gene_objs]
    max_idx = [ idx for idx,value in enumerate(gene_objs_len) if value == max(gene_objs_len)]
    longest_gene = [gene_ids[idx] for idx in max_idx]
    return longest_gene

def genome_pos2trans_pos(data,genome_pos,trans_id):
    '''return the 1_based_position of a site on its transcript based on the genome coordinate'''
    trans_pos = 0
    trans = data.transcript_by_id(trans_id)
    exons = trans.exon_intervals
    strand = trans.strand
    if strand == "+":
        exons = sorted(exons,key=lambda x:x[0])
        for exon in exons:
            if exon[0] <= genome_pos <= exon[1]:
                dis = genome_pos - exon[0] + 1
                trans_pos += dis
            elif genome_pos > exon[1]:
                dis = exon[1] - exon[0] + 1
                trans_pos += dis
            elif exon[0] > genome_pos:
                break
    elif strand == "-":
        exons = sorted(exons,key=lambda x:x[0],reverse=True)
        for exon in exons:
            if exon[0] <= genome_pos <= exon[1]:
                dis = exon[1] - genome_pos + 1
                trans_pos += dis
            elif genome_pos < exon[0]:
                dis = exon[1] - exon[0] +1
                trans_pos += dis
            elif genome_pos > exon[1]:
                break
    return trans_pos

def trans_pos2genome_pos(data,trans_pos,trans_id):
    '''return the genome coordinate of a site base on its transcript coordinate'''
    trans = data.transcript_by_id(trans_id)
    exons = trans.exon_intervals
    genome_pos = 0
    strand = trans.strand
    if strand == "+":
        exons = sorted(exons,key=lambda x:x[0])
        dis = []
        tmp_dis = 0
        for i in range(len(exons)):
            tmp_dis += exons[i][1] - exons[i][0] +1
            dis.append(tmp_dis)
        n_exon = min([idx for idx,value in enumerate(dis) if value >= trans_pos])
        if n_exon > 0:
            genome_pos += exons[n_exon][0] + (trans_pos - dis[n_exon-1] - 1)
        elif n_exon == 0:
            genome_pos += exons[0][0] + (trans_pos - 1)
    elif strand == "-":
        exons = sorted(exons,key=lambda x:x[0],reverse=True)
        dis = []
        tmp_dis = 0
        for i in range(len(exons)):
            tmp_dis += exons[i][1] - exons[i][0] +1
            dis.append(tmp_dis)
        n_exon = min([idx for idx,value in enumerate(dis) if value >= trans_pos])
        if n_exon >0:
            genome_pos += exons[n_exon][1] - (trans_pos - dis[n_exon-1] -1)
        else:
            genome_pos += exons[0][1] - (trans_pos - 1)
    return genome_pos

def get_trans_region(data,trans_id,pos):
    '''return the transcript region of a site based on its trans_pos and its transcript'''
    out = "exon"
    trans = data.transcript_by_id(trans_id)
    strand = trans.strand
    try:
        cds_exons = trans.coding_sequence_position_ranges
        if strand == "+":
            cds_exons = sorted(cds_exons,key=lambda x:x[0])
            if pos < cds_exons[0][0]:
                out = "five_utr"
            elif cds_exons[0][0] <= pos <= cds_exons[-1][1]:
                out = "cds"
            else:
                out = "three_utr"

        elif strand == "-":
            cds_exons = sorted(cds_exons,key=lambda x:x[0],reverse=True)
            if pos > cds_exons[0][1]:
                out = "five_utr"
            elif cds_exons[-1][0] <= pos <= cds_exons[0][1]:
                out = "cds"
            else:
                out = "three_utr"
    except:
        out = "exon"

    return out

def get_five_utr_len(data,trans_id):
    '''return the length of five utr of input coding-transcript id(ENST***)\n
       if input transcript is not protein-coding,return None'''
    trans = data.transcript_by_id(trans_id)
    strand = trans.strand
    exons = trans.exon_intervals
    five_utr_len = 0
    try:
        coding_ranges = trans.coding_sequence_position_ranges
        if strand == "+":
            coding_ranges = sorted(coding_ranges,key=lambda x:x[0])
            exons = sorted(exons,key=lambda x:x[0])
            cds_start = coding_ranges[0][0]
            for i in exons:
                if i[1] < cds_start:
                    five_utr_len += i[1] - i[0] + 1
                elif i[0] < cds_start <= i[1]:
                    five_utr_len += cds_start - i[0]
                elif i[0] >= cds_start:
                    break
        elif strand == "-":
            coding_ranges = sorted(coding_ranges,key=lambda x:x[0],reverse=True)
            exons = sorted(exons,key=lambda x:x[0],reverse=True)
            cds_start = coding_ranges[0][1]
            for i in exons:
                if i[0] > cds_start:
                    five_utr_len += i[1] - i[0] + 1
                elif i[0] <= cds_start < i[1]:
                    five_utr_len += i[1] - cds_start
                elif i[1] <= cds_start:
                    break
        return five_utr_len
    except:
        print("%s is not coding transcript!" % trans_id)
        return None
