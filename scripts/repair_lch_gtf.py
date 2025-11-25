from collections import namedtuple, OrderedDict
import sys

# Parse arguments
args = sys.argv[1:]
input_gff, output_gff = args

# Utilities
GTF_record = namedtuple("gtf_record", [
    "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"
])

def gtf_line_parser(gtf_string):
    """Parse a gtf line into a GTF_record namedtuple."""
    fields = gtf_string.strip().split("\t")
    attrs = {}
    for attr in fields[8].split(";"):
        if attr.strip() == "":
            continue
        key, value = attr.strip().split(" ", 1) # split on the first space
        attrs[key] = value.strip().strip('"')
    
    return GTF_record(
        seqname=fields[0],
        source=fields[1],
        feature=fields[2],
        start=int(fields[3]),
        end=int(fields[4]),
        score=fields[5],
        strand=fields[6],
        frame=fields[7],
        attribute=attrs
    )
    
def gtf_tuple_to_line(gtf_record):
    txt = "\t".join([str(i) for i in gtf_record[:8]]) + "\t"
    txt += ";".join([f'{k} "{v}"' for k, v in gtf_record.attribute.items()]) + ";\n"
    return txt
    
def gtf_to_dict(fname):
    gtf_dict = OrderedDict()
    with open(fname, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            record = gtf_line_parser(line)
            transcript_entries = gtf_dict.get(record.attribute["transcript_id"], [])
            transcript_entries.append(record)
            gtf_dict[record.attribute["transcript_id"]] = transcript_entries
    
    return gtf_dict

# Import GTF as dict 
gtf_dict = gtf_to_dict(input_gff)

# Repair GTF
for k, v in gtf_dict.items():
    # add exon entries if missing
    has_exon = any([entry.feature == "exon" for entry in v])
    if not has_exon:
        # Find all CDS entries
        cds_entries = [entry for entry in v if entry.feature == "CDS"]
        # For each CDS entry, create a corresponding exon entry
        for cds in cds_entries:
            exon_entry = GTF_record(
                seqname=cds.seqname,
                source="repaire",
                feature="exon",
                start=cds.start,
                end=cds.end,
                score=cds.score,
                strand=cds.strand,
                frame=".",
                attribute=cds.attribute
            )
            gtf_dict[k].append(exon_entry)
        
        # Sort entries by start position
        gtf_dict[k] = sorted(gtf_dict[k], key=lambda x: x.start)
    
    # add gene_id attribute if missing
    for entry in v:
        if "gene_id" not in entry.attribute:
            entry.attribute["gene_id"] = entry.attribute["transcript_id"]


with open(output_gff, "w") as out_f:
    for tid, records in gtf_dict.items():
        for record in records:
            foo = out_f.write(gtf_tuple_to_line(record))
