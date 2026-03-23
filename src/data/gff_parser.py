import gzip
from pathlib import Path
from Bio import SeqIO
from io import StringIO


def open_gff_as_text(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "r")

def extract_gene_seq(gff_file_path, gene_name:str )->str:
    with open_gff_as_text(gff_file_path) as f:
        content = f.read()

        if "##FASTA" not in content:
            return None

        annotation_section, fasta_section = content.split("##FASTA")
        gene_coords = None

        for line in annotation_section.splitlines():
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = cols
            if ftype != "CDS":
                continue

            attrs_dict = {}
            for attr in attrs.split(";"):
                if "=" in attr:
                    key,val = attr.split("=")
                    attrs_dict[key.strip().lower()] = val.strip().lower()
            gene_id = attrs_dict.get("name", "")
            locus_tag = attrs_dict.get("locus_tag", "")
            product = attrs_dict.get("product", "")
            if (gene_name.lower() in gene_id.lower() or
                    gene_name.lower() in locus_tag.lower() or
                    gene_name.lower() in product.lower()):
                gene_coords = {
                    "seqid": seqid,
                    "start": int(start),  # Convert string to int
                    "end": int(end),
                    "strand": strand
                }
                break  # Stop searching once found

        if not gene_coords:
            return None

        fasta_handle = StringIO(fasta_section.lstrip())
        for record in SeqIO.parse(fasta_handle, "fasta"):
            if record.id == gene_coords["seqid"]:
                start_idx = gene_coords["start"] - 1
                end_idx = gene_coords["end"]
                gene_seq = record.seq[start_idx:end_idx]
                if gene_coords["strand"] == "-":
                    gene_seq = gene_seq.reverse_complement()
                return str(gene_seq)
        return None


def extract_fasta(gff_file_path):
    with open_gff_as_text(gff_file_path) as f:
        content = f.read()

        if "##FASTA" not in content:
            return None

        annotation_section, fasta_section = content.split("##FASTA")
        fasta_handle = StringIO(fasta_section.lstrip())
        return SeqIO.parse(fasta_handle, "fasta")


def parse_attrs(attr_field: str) -> dict:
    d = {}
    for part in attr_field.split(";"):
        if not part:
            continue
        k, _, v = part.partition("=")
        d[k] = v
    return d


def scan_ko_hits(gff_path: Path, target_kos: set) -> set:
    ko_hits = set()
    with open_gff_as_text(gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype, attrs = cols[2], cols[8]
            if ftype not in ("CDS", "gene"):
                continue
            d = parse_attrs(attrs)
            kegg_field = d.get("KEGG", "")
            if kegg_field and kegg_field != "-":
                ko_entries = [ko.strip() for ko in kegg_field.split(",")]
                for ko in ko_entries:
                    ko_id = ko.replace("ko:","")
                    if ko_id in target_kos:
                        ko_hits.add(ko_id)
    return ko_hits

def extract_eggnog_id(gff_path: Path):
    eggnog_ids = set()
    with open_gff_as_text(gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype, attrs = cols[2], cols[8]
            if ftype not in ("CDS", "gene"):
                continue
            d = parse_attrs(attrs)
            eggnog_field = d.get("eggNOG", "")
            if eggnog_field and eggnog_field != "-":
                for eggnog_id in eggnog_field.split(","):
                    eggnog_ids.add(eggnog_id.strip())
    return eggnog_ids