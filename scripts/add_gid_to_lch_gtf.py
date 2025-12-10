import sys

args = sys.argv[1:]
input_gtf, output_gtf = args

with open(input_gtf, "r") as f, open(output_gtf, "w") as g:
    for line in f:
        if line.startswith("#"):
            g.write(line)
        else:
            # add gene_id {id} if missing gene_id
            fields = line.split("\t")
            if not "gene_id" in fields[-1]:
                extended_attribute = fields[-1].replace(
                    "transcript_id", "gene_id"
                )
                fields[-1] = fields[-1].strip() + " " + extended_attribute
                
            g.write("\t".join(fields))