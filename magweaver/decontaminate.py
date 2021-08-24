import os
from collections import defaultdict, Counter
from magweaver import genome

BASEPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP_DIR = os.path.join(BASEPATH, "tmp")
OUT_DIR = os.path.join(BASEPATH + "results")

def create_mag(mag_fasta, forward_reads, reverse_reads, tmp_dir = TMP_DIR):
    mag = genome.Mag(mag_fasta, forward_reads, reverse_reads, tmp_dir)
    mag.craft_mag()
    mag.craft_summary_dict()
    return mag

def flag_erroneous_contigs(mag, minimum_suspicion):
    all_lists = flag_gc(mag) + flag_cov(mag) + flag_tetra(mag) + flag_tax(mag) + flag_scg_duplicates(mag) + flag_scg_empty(mag) + flag_mobilome(mag)
    suspicion_dict = dict(Counter(all_lists))
    highly_suspect_contigs = []
    for contig, suspicion in suspicion_dict.items():
        if suspicion >= minimum_suspicion:
            highly_suspect_contigs.append(contig)
    return highly_suspect_contigs           

def filter_mag(mag, highly_suspect_contigs):
    filtered_mag = {contig: seqcontents for contig, seqcontents in mag.mag_contigs.items() if contig not in highly_suspect_contigs}
    return filtered_mag

def write_fasta(filtered_mag, outfile):
    for contig, seqcontents in filtered_mag.items():
        contig_name = contig
        contig_seq = str(seqcontents["seq"])
        with open(outfile, 'a+') as outf:
            outf.write(">" + str(contig_name) + "\n" + contig_seq + "\n")

    
# ------ Feeder functions for flag_erroneous_contigs ------   
def flag_gc(mag):
    erroneous_contigs = []
    min_gc = mag.summary_dict["gc_mean"] - mag.summary_dict["gc_std"]
    max_gc = mag.summary_dict["gc_mean"] - mag.summary_dict["gc_std"]
    for contig, seqcontents in mag.mag_contigs.items():
        gc = seqcontents["gc"]
        if (gc > max_gc) or (gc < min_gc):
            erroneous_contigs.append(contig)
    return erroneous_contigs

def flag_cov(mag):
    erroneous_contigs = []
    min_cov = mag.summary_dict["cov_mean"] - mag.summary_dict["cov_std"]
    max_cov = mag.summary_dict["cov_mean"] + mag.summary_dict["cov_std"]
    for contig, seqcontents in mag.mag_contigs.items():
        cov = seqcontents["cov"]
        if (cov > max_cov) or (cov < min_cov):
            erroneous_contigs.append(contig)
    return erroneous_contigs

def flag_tetra(mag):
    erroneous_contigs = []
    min_comp = mag.summary_dict["pca_mean"] - mag.summary_dict["pca_std"]
    max_comp = mag.summary_dict["pca_mean"] + mag.summary_dict["pca_std"]
    for contig, seqcontents in mag.mag_contigs.items():
        comp = seqcontents["pca"]
        if (comp > min_comp) or (comp < max_comp):
            erroneous_contigs.append(contig)
    return erroneous_contigs

def flag_tax(mag):
    erroneous_contigs = []
    consensus_phylum = mag.summary_dict["phylum"]
    for contig, seqcontents in mag.mag_contigs.items():
        if seqcontents["tax"]["phylum"] != consensus_phylum:
            erroneous_contigs.append(contig)
    return erroneous_contigs

def flag_scg_duplicates(mag):
    potential_duplicates = []
    scg_dict = defaultdict(int)
    for seqcontents in mag.mag_contigs.values():
        contig_scg_dict = seqcontents["scg"]
        for hmm in contig_scg_dict.keys():
            scg_dict[hmm] += 1
    duplicates = {hmm: hits for hmm, hits in scg_dict.items() if hits > 1}
    for contig, seqcontents in mag.mag_contigs.items():
        for hmm in duplicates.keys():
            if hmm in seqcontents["scg"].keys():
                potential_duplicates.append(contig)
    return potential_duplicates

def flag_scg_empty(mag):
    no_hits = []
    scg_dict = defaultdict(int)
    for contig, seqcontents in mag.mag_contigs.items():
        if not seqcontents["scg"]:
            no_hits.append(contig)
    return no_hits
    
def flag_mobilome(mag):
    erroneous_contigs = []
    for contig, seqcontents in mag.mag_contigs.items():
        if seqcontents["mob"]:
            erroneous_contigs.append(contig)
    return erroneous_contigs        

# ------ IT'S TIME FOR THE DECONTAMINATOR ------

def decontaminator(mag_fasta, forward_reads, reverse_reads, minimum_suspicion, outfile, tmp_dir = TMP_DIR):
    if not os.path.exists(OUT_DIR):
        os.mkdir(OUT_DIR)
    outloc = os.path.join(OUT_DIR, outfile)
    mag = create_mag(mag_fasta, forward_reads, reverse_reads, tmp_dir)
    suspect_contigs = flag_erroneous_contigs(mag, minimum_suspicion)
    filtered_mag = filter_mag(mag, suspect_contigs)
    write_fasta(filtered_mag, outloc)
    