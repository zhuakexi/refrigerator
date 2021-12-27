import bioframe
def get_arms(ref_name,clr):
    chromsizes = bioframe.fetch_chromsizes(ref_name)
    cens = bioframe.fetch_centromeres(ref_name)
    arms = bioframe.core.construction.add_ucsc_name_column(bioframe.make_chromarms(chromsizes, cens))
    arms = arms[arms.chrom.isin(clr.chromnames)].reset_index(drop=True)
    return arms