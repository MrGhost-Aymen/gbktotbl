#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, BeforePosition, AfterPosition
import sys
import os
import logging
from typing import Dict, List, Tuple, Optional, Set, DefaultDict
from collections import defaultdict
import argparse

# Configure logging
logging.basicConfig(level=logging.WARNING, format='%(levelname)s: %(message)s')

__version__ = "1.5.0"

# Valid INSDC/NCBI features
VALID_FEATURES = {
    "gene", "CDS", "mRNA", "tRNA", "rRNA", "ncRNA", "exon", "intron",
    "5'UTR", "3'UTR", "misc_feature", "repeat_region", "D-loop"
}

# Features that do not require CDS counterparts
NON_CDS_FEATURES = {"tRNA", "rRNA", "ncRNA", "misc_feature", "repeat_region", "D-loop"}

# Protein-coding gene prefixes (used to detect likely coding genes)
PROTEIN_CODING_GENE_PREFIXES = {"nad", "cox", "atp", "cytb", "mat", "sdh", "ccm", "rps", "rpl"}

# Recommended qualifier order (NCBI best practice)
QUALIFIER_ORDER = [
    'gene', 'locus_tag', 'function', 'product', 'protein_id',
    'transcript_id', 'note', 'codon_start', 'transl_table'
]

# Qualifiers that require special handling
QUALIFIER_SANITIZATION = {
    'gene': lambda x: x.replace(" ", "_").replace("(", "").replace(")", ""),
    'product': lambda x: x if not x.isupper() else x.capitalize(),
    'note': lambda x: x.replace("\n", " ").replace("\t", " ").strip(),
}


def setup_logging(verbose: bool = False) -> None:
    """Configure logging level based on verbosity."""
    level = logging.DEBUG if verbose else logging.WARNING
    logging.basicConfig(level=level, format='%(levelname)s: %(message)s')


def sanitize_qualifier_value(key: str, value: str) -> str:
    """Sanitize qualifier values according to NCBI standards."""
    if key in QUALIFIER_SANITIZATION:
        return QUALIFIER_SANITIZATION[key](value)
    return value.strip()


def wrap_qualifier(text: str, width: int = 60) -> str:
    """Wrap qualifier text at specified width for TBL output."""
    lines = []
    for i in range(0, len(text), width):
        lines.append(text[i:i+width])
    return "\n\t\t\t\t".join(lines)


def is_partial(location: FeatureLocation) -> bool:
    """Check if a feature location is partial (5' or 3')."""
    if isinstance(location, CompoundLocation):
        return (isinstance(location.parts[0].start, BeforePosition) or
                isinstance(location.parts[-1].end, AfterPosition))
    return (isinstance(location.start, BeforePosition) or
            isinstance(location.end, AfterPosition))


def calculate_codon_start(cds_seq: str) -> Optional[str]:
    """Calculate proper codon_start for partial CDS features."""
    seq_len = len(cds_seq)
    if seq_len % 3 == 0:
        return None
    frame = (3 - seq_len % 3) % 3
    return str(frame + 1)  # Returns 1, 2, or 3


def get_codon_start(feature, sequence) -> Optional[str]:
    """Infer codon_start for partial CDS based on sequence length."""
    if feature.type != "CDS" or not is_partial(feature.location):
        return None

    try:
        cds_seq = feature.location.extract(sequence).seq
        return calculate_codon_start(str(cds_seq))
    except Exception as e:
        logging.warning(f"Could not compute codon_start for CDS in {feature}: {e}")
        return None


def format_location(feature) -> List[Tuple[str, str]]:
    """Convert BioPython location to NCBI TBL format with proper sorting."""
    location = feature.location

    if isinstance(location, CompoundLocation):
        # Sort parts by start position
        location.parts.sort(key=lambda loc: loc.start)
        parts = location.parts
    else:
        parts = [location]

    formatted = []
    for part in parts:
        start = part.start + 1  # Convert to 1-based
        end = part.end

        start_str = f"<{start}" if isinstance(part.start, BeforePosition) else str(start)
        end_str = f">{end}" if isinstance(part.end, AfterPosition) else str(end)

        # Handle reverse strand
        if part.strand == -1:
            start_str, end_str = end_str, start_str

        formatted.append((start_str, end_str))

    return formatted


def validate_sequence(record) -> None:
    """Check sequence for potential submission issues."""
    seq = str(record.seq)
    if "N" in seq:
        n_count = seq.count("N")
        logging.warning(f"Record {record.id} contains {n_count} Ns in sequence")
    if len(seq) < 50:
        logging.warning(f"Record {record.id} is very short ({len(seq)} bp)")


def check_gene_cds_coverage(gene_feature, cds_list):
    """Check if gene covers all intervals of its CDS features."""
    gene_loc = gene_feature.location
    gene_start = gene_loc.start
    gene_end = gene_loc.end

    for cds in cds_list:
        cds_loc = cds.location
        cds_start = cds_loc.start
        cds_end = cds_loc.end

        # For reverse strand, compare positions directly
        if cds_start < gene_start or cds_end > gene_end:
            logging.warning(
                f"Gene-CDS location mismatch for '{gene_feature.qualifiers.get('gene', [''])[0]}': "
                f"gene({gene_start + 1}-{gene_end}) vs CDS({cds_start + 1}-{cds_end}). "
                f"Gene should span entire CDS range."
            )


def check_feature_consistency(record) -> None:
    """Check consistency between gene and CDS/mRNA/tRNA features."""
    genes = {}
    cdses_by_gene = defaultdict(list)
    features_by_gene = defaultdict(list)

    for feature in record.features:
        if feature.type == "gene":
            gene_id = feature.qualifiers.get('gene', [''])[0]
            if gene_id:
                genes[gene_id] = feature
        elif feature.type in ["CDS", "mRNA"]:
            gene_ids = feature.qualifiers.get('gene', [])
            for gene_id in gene_ids:
                if feature.type == "CDS":
                    cdses_by_gene[gene_id].append(feature)
                features_by_gene[gene_id].append(feature)
        elif feature.type in NON_CDS_FEATURES:
            gene_ids = feature.qualifiers.get('gene', [])
            for gene_id in gene_ids:
                features_by_gene[gene_id].append(feature)

    # Now validate each gene
    for gene_id, gene_feature in genes.items():
        if gene_id in NON_CDS_FEATURES:
            continue  # Skip known non-CDS features

        if gene_id not in features_by_gene:
            logging.warning(f"Gene '{gene_id}' has no associated features (CDS, mRNA, etc.)")

        if gene_id not in cdses_by_gene:
            # Only warn if it's likely a protein-coding gene
            if any(gene_id.lower().startswith(prefix) for prefix in PROTEIN_CODING_GENE_PREFIXES):
                logging.warning(f"Gene '{gene_id}' appears protein-coding but has no CDS")
            continue

        # Check if gene fully covers all CDS intervals
        check_gene_cds_coverage(gene_feature, cdses_by_gene[gene_id])


def write_feature_qualifiers(tbl_out, feature, suppress_gene: bool = False) -> None:
    """Write properly formatted qualifiers to TBL file in recommended order."""
    written_qualifiers = set()

    # Suppress overlapping gene qualifier
    gene_name = feature.qualifiers.get('gene', [''])[0]
    if gene_name and suppress_gene:
        tbl_out.write(f"\t\t\tgene\t-\n")
        written_qualifiers.add('gene')

    # First write qualifiers in recommended order
    for key in QUALIFIER_ORDER:
        if key in feature.qualifiers:
            written_qualifiers.add(key)
            for value in feature.qualifiers[key]:
                for v in value.split(","):
                    v = sanitize_qualifier_value(key, v.strip())
                    if v:
                        wrapped = wrap_qualifier(v)
                        tbl_out.write(f"\t\t\t{key}\t{wrapped}\n")

    # Then write remaining qualifiers not in the recommended order
    for key in sorted(feature.qualifiers.keys()):
        if key not in written_qualifiers and key != "translation":
            for value in feature.qualifiers[key]:
                for v in value.split(","):
                    v = sanitize_qualifier_value(key, v.strip())
                    if v:
                        wrapped = wrap_qualifier(v)
                        tbl_out.write(f"\t\t\t{key}\t{wrapped}\n")


def write_tbl(genbank_file: str, tbl_file: str,
              transl_table: Optional[int] = None,
              reference_pmid: Optional[str] = None,
              check_consistency: bool = True,
              offset: Optional[int] = None) -> None:
    """Convert GenBank to NCBI TBL format with strict validation."""
    global_offset = offset if offset is not None else 0

    try:
        if not os.path.isfile(genbank_file):
            raise FileNotFoundError(f"Input file '{genbank_file}' does not exist")

        with open(tbl_file, "w") as tbl_out:
            for record in SeqIO.parse(genbank_file, "genbank"):
                validate_sequence(record)
                if check_consistency:
                    check_feature_consistency(record)

                # Write header
                tbl_out.write(f">Feature {record.id}\n")

                # Add REFERENCE line
                if reference_pmid:
                    tbl_out.write(f"1\t{len(record)}\tREFERENCE\n")
                    tbl_out.write(f"\t\t\tPubMed\t{reference_pmid}\n")

                # Add offset
                if offset:
                    tbl_out.write(f"[offset={offset}]\n")

                for feature in record.features:
                    if feature.type == "source":
                        continue

                    if feature.type not in VALID_FEATURES:
                        logging.warning(f"Skipping non-standard feature '{feature.type}' in {record.id}")
                        continue

                    # Validate CDS features
                    if feature.type == "CDS":
                        if "product" not in feature.qualifiers:
                            logging.warning(f"Skipping CDS in {record.id} missing 'product' qualifier")
                            continue
                        # Optionally add transl_table
                        if transl_table:
                            feature.qualifiers.setdefault("transl_table", [str(transl_table)])

                    # Process feature locations
                    location_parts = format_location(feature)

                    # Write feature lines
                    for i, (start, end) in enumerate(location_parts):
                        tbl_out.write(f"{start}\t{end}\t{feature.type}\n")

                        if i == len(location_parts) - 1:  # Only for last part
                            # Add codon_start if needed
                            codon_start = get_codon_start(feature, record)
                            if codon_start:
                                tbl_out.write(f"\t\t\tcodon_start\t{codon_start}\n")

                            # Write qualifiers
                            write_feature_qualifiers(tbl_out, feature, suppress_gene=True)

    except Exception as e:
        logging.error(f"Fatal error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert GenBank files to NCBI TBL format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gbk_to_tbl.py input.gbk output.tbl --transl_table 2 --pmid 12345678
  python gbk_to_tbl.py input.gbk output.tbl --no-consistency-check -v
        """
    )
    parser.add_argument("input", help="Input GenBank file (.gbk)")
    parser.add_argument("output", help="Output TBL file")
    parser.add_argument("--transl_table", type=int, default=None,
                        help="Add transl_table (e.g., 2 for mitochondria)")
    parser.add_argument("--pmid", dest="pmid", default=None,
                        help="Add PubMed ID as a REFERENCE feature")
    parser.add_argument("--no-consistency-check", action="store_true",
                        help="Skip gene/CDS consistency checks")
    parser.add_argument("--offset", type=int, default=None,
                        help="Add an offset to all positions after this file")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable debug-level logging")
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {__version__}")

    args = parser.parse_args()

    setup_logging(args.verbose)

    try:
        # Test output file accessibility
        with open(args.output, "w") as f:
            pass
    except IOError as e:
        logging.error(f"Cannot write to '{args.output}': {e}")
        sys.exit(1)

    write_tbl(args.input, args.output,
              transl_table=args.transl_table,
              reference_pmid=args.pmid,
              check_consistency=not args.no_consistency_check,
              offset=args.offset)

    logging.info(f"Successfully converted {args.input} to {args.output}")