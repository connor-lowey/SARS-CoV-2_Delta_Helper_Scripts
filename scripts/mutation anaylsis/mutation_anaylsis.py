import csv
import math
from collections import Counter

import pandas as pd

CONFIG_PATH = "mutation_analysis_config.csv"
SIGNIFICANCE_THRESHOLD_PCT = 1

REGION_COLUMNS = ["Position", "A", "G", "C", "U", "GAP"]
MUTATION_COLUMNS = ["A", "G", "C", "U", "GAP"]
BASE_TO_COLUMN = {"A": "A", "G": "G", "C": "C", "T": "U", "U": "U", "-": "GAP"}
BASE_TO_RNA = {"A": "A", "G": "G", "C": "C", "T": "U", "U": "U", "-": "-"}
COLUMN_TO_BASE = {"A": "A", "G": "G", "C": "C", "U": "U", "GAP": "-"}

CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "---": "-",
}


def read_config(config_path=CONFIG_PATH):
    """Read the config CSV, returning (data_file, msa_file, genomes, group, regions).

    The file starts with "label,value" metadata rows (positions file,
    alignment file, genome count, group name), followed by a blank line,
    then "name,start,end" region rows.
    """
    with open(config_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        rows = [row for row in reader]

    blank_index = rows.index([])
    metadata_rows = [row for row in rows[:blank_index] if row]
    region_rows = rows[blank_index + 1:]

    # metadata rows are positional: positions file, alignment file, genome count, group
    data_file = metadata_rows[0][1]
    msa_file = metadata_rows[1][1]
    genomes = int(metadata_rows[2][1])
    group = metadata_rows[3][1]

    regions = {}
    for row in region_rows:
        if not row:
            continue
        name, start, end = row[0], row[1], row[2]
        regions[name] = (int(start), int(end))

    return data_file, msa_file, genomes, group, regions


def read_position_summary(data_file):
    """Read the tab-delimited nt-position summary txt file into a DataFrame.

    The file has a couple of preamble lines (title, alignment path, blank
    line) before the "Position A G C U GAP" header, so locate that header
    line rather than assuming a fixed number of rows to skip.
    """
    with open(data_file, encoding="utf-8-sig") as f:
        lines = f.readlines()

    header_index = next(i for i, line in enumerate(lines) if line.strip().startswith("Position"))
    return pd.read_csv(data_file, sep="\t", skiprows=header_index, encoding="utf-8-sig")


def split_by_region(data_file, regions):
    """Return a dict mapping sheet name to DataFrame: the whole genome
    plus one entry per region, kept in memory instead of written to xlsx.
    """
    whole_genome = read_position_summary(data_file)
    position_col = whole_genome.columns[0]

    split_data = {"Whole_genome": whole_genome}

    for name, (start, end) in regions.items():
        mask = (whole_genome[position_col] >= start) & (whole_genome[position_col] <= end)
        region_df = whole_genome.loc[mask].copy()
        region_df.columns = REGION_COLUMNS[: len(region_df.columns)]
        split_data[name] = region_df

    return split_data


def read_reference_sequence(fasta_path):
    """Return the first (reference) sequence from an MSA fasta file, gaps included."""
    return read_all_sequences(fasta_path)[0]


def read_all_sequences(fasta_path):
    """Return every aligned sequence in an MSA fasta file, reference first, gaps included."""
    with open(fasta_path) as f:
        lines = f.read().splitlines()

    sequences = []
    current = []
    for line in lines:
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
            current = []
        else:
            current.append(line)
    if current:
        sequences.append("".join(current))

    return sequences


def write_split_workbook(split_data, output_file):
    """Write the whole genome sheet plus one sheet per region to an xlsx file."""
    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
        for name, df in split_data.items():
            df.to_excel(writer, sheet_name=name[:31], index=False)


def build_region_frame(reference_seq, start, end):
    """Map alignment columns in [start, end] to their 1-based rank among
    non-gap reference bases (the gene position with reference gaps ignored).

    Used to report the gap-excluded position columns, and to compute the
    codon for non-insertion mutations against the true, ungapped reading
    frame -- an upstream reference-gap run that isn't a multiple of 3
    otherwise shifts raw alignment-column arithmetic out of frame. Gap
    columns map to the rank of the next real base (best-effort, since an
    inserted base has no ungapped position of its own).
    """
    column_to_real = {}
    real_to_column = []
    for col in range(start, end + 1):
        if reference_seq[col - 1] != "-":
            real_to_column.append(col)
            column_to_real[col] = len(real_to_column)
        else:
            column_to_real[col] = len(real_to_column) + 1
    return column_to_real, real_to_column


def get_codon_context(reference_seq, real_to_column, real_position):
    """Return (offset, codon_positions, original_codon) for the codon containing
    real_position, using the ungapped reference frame so an upstream reference
    gap run that isn't a multiple of 3 doesn't shift the codon out of frame.
    """
    offset = (real_position - 1) % 3
    codon_start_index = real_position - offset
    codon_positions = real_to_column[codon_start_index - 1: codon_start_index + 2]
    original_codon = [BASE_TO_RNA.get(reference_seq[pos - 1].upper(), reference_seq[pos - 1].upper()) for pos in codon_positions]
    return offset, codon_positions, original_codon


def get_insertion_codon_context(reference_seq, start, before_nt, before_column, full_position):
    """Return (offset, codon_positions, original_codon) for an insertion's own
    local codon. The reading frame continues from the real reference codon
    immediately preceding the insertion (before_nt/before_column), so an
    upstream reference-gap run that isn't a multiple of 3 shifts this frame
    exactly like it does for real (non-insertion) positions. Falls back to
    naive local arithmetic if there's no real base before the insertion.
    """
    if before_column is None:
        offset = (full_position - start) % 3
    else:
        offset = (before_nt + (full_position - before_column - 1)) % 3
    codon_start = full_position - offset
    codon_positions = [codon_start, codon_start + 1, codon_start + 2]
    original_codon = [BASE_TO_RNA.get(reference_seq[pos - 1].upper(), reference_seq[pos - 1].upper()) for pos in codon_positions]
    return offset, codon_positions, original_codon


def get_insertion_flanks(real_to_column, full_position, column_to_real):
    """Return (before_nt, after_nt, before_column): the ungapped gene positions
    flanking an insertion, and the full genome column of the real reference
    base immediately preceding it.
    """
    after_nt = column_to_real[full_position]
    before_nt = after_nt - 1
    before_column = real_to_column[before_nt - 1] if before_nt >= 1 else None
    return before_nt, after_nt, before_column


def format_flank(before, after):
    """Format a flanking pair as "before/after", collapsing to a single
    value when both sides land on the same position (e.g. same codon).
    """
    if before == after:
        return str(before)
    return f"{before}/{after}"


def build_significant_bases(region_df):
    """Map each position to the set of mutant base letters occurring in at
    least 1% of genomes there (the same threshold used to decide whether a
    mutation gets its own reported row).

    Used to normalize noise when building the majority codon: a nearby
    variant that never clears the reporting threshold on its own shouldn't
    be treated as a genuine co-occurring mutation. Uses each position's own
    included count as the denominator, matching the main reporting filter.
    """
    significant = {}
    for _, row in region_df.iterrows():
        position = int(row["Position"])
        included = sum(int(row[col]) for col in MUTATION_COLUMNS) - 1
        bases = set()
        for col in MUTATION_COLUMNS:
            if int(row[col]) / included * 100 >= SIGNIFICANCE_THRESHOLD_PCT:
                bases.add(COLUMN_TO_BASE[col])
        significant[position] = bases
    return significant


def get_codon_change(all_sequences, offset, codon_positions, original_codon, mutant_column, significant_bases):
    """Return "<original codon>-<mutated codon>" for the codon containing full_position.

    The mutated codon is the most common actual codon among genomes that
    carry this mutation -- i.e. whichever is more prevalent between carrying
    it alone (other codon positions at reference) vs. together with another
    mutation elsewhere in the same codon. At each other position, any base
    that doesn't independently clear the 1% reporting threshold there is
    normalized to the reference base first, so a rare/noise variant at a
    neighboring position doesn't get mistaken for a genuine co-occurring
    mutation.
    """
    mutant_base = COLUMN_TO_BASE[mutant_column]

    codon_counts = Counter()
    for seq in all_sequences:
        raw_bases = [BASE_TO_RNA.get(seq[pos - 1].upper(), seq[pos - 1].upper()) for pos in codon_positions]
        if raw_bases[offset] != mutant_base:
            continue

        normalized = tuple(
            base if base in significant_bases.get(pos, ()) else original_codon[i]
            for i, (pos, base) in enumerate(zip(codon_positions, raw_bases))
        )
        codon_counts[normalized] += 1

    if codon_counts:
        mutated_codon = list(codon_counts.most_common(1)[0][0])
    else:
        mutated_codon = list(original_codon)
        mutated_codon[offset] = mutant_base

    return f"{''.join(original_codon)}-{''.join(mutated_codon)}"


def get_codon_cooccurrence(all_sequences, offset, codon_positions, original_codon, mutant_column, included, significant_bases):
    """Return a summary string splitting genomes with this mutation into those
    that carry it alone vs. together with another mutation elsewhere in the
    same codon (in the same genome). Percentages are out of `included`, the
    same per-position denominator used for the main percentage column.
    """
    mutant_base = COLUMN_TO_BASE[mutant_column]
    alone = 0
    together = 0

    for seq in all_sequences:
        base_at_offset = BASE_TO_RNA.get(seq[codon_positions[offset] - 1].upper(), seq[codon_positions[offset] - 1].upper())
        if base_at_offset != mutant_base:
            continue

        has_other_mutation = False
        for other_offset in range(3):
            if other_offset == offset:
                continue
            other_pos = codon_positions[other_offset]
            other_base = BASE_TO_RNA.get(seq[other_pos - 1].upper(), seq[other_pos - 1].upper())
            if other_base not in significant_bases.get(other_pos, ()):
                continue
            if other_base != original_codon[other_offset]:
                has_other_mutation = True
                break

        if has_other_mutation:
            together += 1
        else:
            alone += 1

    alone_pct = round(alone / included * 100, 2)
    together_pct = round(together / included * 100, 2)
    return f"Alone: {alone} ({alone_pct}%); Together: {together} ({together_pct}%)"


def get_amino_acid_change(codon_change):
    """Return "<original AA>-<mutated AA>" for a "<original codon>-<mutated codon>" string.

    Codons containing a gap or other non-standard base (e.g. from a partial
    deletion) can't be translated, so "?" is used for that side. Sliced by
    fixed position rather than split("-") since a codon can itself contain
    a gap character ("-").
    """
    original_codon, mutated_codon = codon_change[:3], codon_change[4:7]
    original_aa = CODON_TABLE.get(original_codon, "?")
    mutated_aa = CODON_TABLE.get(mutated_codon, "?")
    return f"{original_aa}-{mutated_aa}"


def get_mutation_type(ref_column, mutant_column, amino_acid_change):
    """Classify a mutation as Insertion (reference itself is a gap here),
    Deletion (the mutant nucleotide is a gap), or Synonymous/Non-synonymous
    (based on the amino acid change) otherwise.
    """
    if ref_column == "GAP":
        return "Insertion"
    if mutant_column == "GAP":
        return "Deletion"
    return "Synonymous" if amino_acid_change[0] == amino_acid_change[-1] else "Non-synonymous"


def build_results(split_data, regions, all_sequences, genomes, group):
    """Build the mutation results table.

    One row is produced per (region position, mutant nucleotide) where the
    position has more than one non-zero column among A/G/C/U/GAP.
    """
    reference_seq = all_sequences[0]
    percentage_column = f"Percentage of group affected (occurance / number of genomes included)"
    records = []

    for region_name, (start, end) in regions.items():
        region_df = split_data[region_name]
        column_to_real, real_to_column = build_region_frame(reference_seq, start, end)
        significant_bases = build_significant_bases(region_df)

        for _, row in region_df.iterrows():
            full_position = int(row["Position"])
            counts = {col: int(row[col]) for col in MUTATION_COLUMNS}
            nonzero_columns = [col for col, value in counts.items() if value != 0]

            if len(nonzero_columns) <= 1:
                continue

            included = sum(counts.values()) - 1 # Subtract 1 to exclude the reference genome from the included count
            excluded = genomes - included

            ref_base = reference_seq[full_position - 1].upper()
            ref_column = BASE_TO_COLUMN.get(ref_base, ref_base)

            alignment_gene_position_nt = full_position - start + 1  # gaps counted as real positions

            if ref_column == "GAP":
                before_nt, after_nt, before_column = get_insertion_flanks(real_to_column, full_position, column_to_real)
                gene_position_nt = format_flank(before_nt, after_nt)
                gene_position_aa = format_flank(math.ceil(before_nt / 3), math.ceil(after_nt / 3))
                full_genome_position = format_flank(start - 1 + before_nt, start - 1 + after_nt)
                offset, codon_positions, original_codon = get_insertion_codon_context(reference_seq, start, before_nt, before_column, full_position)
            else:
                real_position = column_to_real[full_position]
                gene_position_nt = real_position
                gene_position_aa = math.ceil(real_position / 3)
                full_genome_position = start - 1 + real_position  # shift full_position back by the gaps seen so far in this region
                offset, codon_positions, original_codon = get_codon_context(reference_seq, real_to_column, real_position)

            for mutant_column in nonzero_columns:
                if mutant_column == ref_column:
                    continue

                occurrence = counts[mutant_column]
                percentage = round(occurrence / included * 100, 2)
                if percentage < SIGNIFICANCE_THRESHOLD_PCT:
                    continue

                codon_change = get_codon_change(all_sequences, offset, codon_positions, original_codon, mutant_column, significant_bases)
                codon_cooccurrence = get_codon_cooccurrence(all_sequences, offset, codon_positions, original_codon, mutant_column, included, significant_bases)
                amino_acid_change = get_amino_acid_change(codon_change)

                records.append({
                    "Group": group,
                    "Genome region": region_name,
                    "Gene position (nucleotide) - Reference based": gene_position_nt,
                    "Gene position (amino acid) - Reference based": gene_position_aa,
                    "Full genome position (nucleotide) - Reference based": full_genome_position,
                    "Number of genomes included": included,
                    "Number of genomes excluded": excluded,
                    "Mutation occurrence count": occurrence,
                    percentage_column: f"{percentage:.2f}",
                    "Original nucleotide (in reference genome)": ref_column,
                    "Mutated nucleotide": mutant_column,
                    "Codon Change": codon_change,
                    "Amino Acid Change": amino_acid_change,
                    "Mutation Type": get_mutation_type(ref_column, mutant_column, amino_acid_change),
                    "Codon co-occurrence (alone vs together)": codon_cooccurrence,
                    "Gene position (nucleotide) - Alignment based": alignment_gene_position_nt,
                    "Full genome position (nucleotide) - Alignment based": full_position,
                })

    # Only meaningful when another row shares the same codon (amino acid position).
    aa_position_counts = Counter((record["Genome region"], record["Gene position (amino acid) - Reference based"]) for record in records)
    for record in records:
        key = (record["Genome region"], record["Gene position (amino acid) - Reference based"])
        if aa_position_counts[key] <= 1:
            record["Codon co-occurrence (alone vs together)"] = ""

    return pd.DataFrame(records)


if __name__ == "__main__":
    data_file, msa_file, genomes, group, regions = read_config()
    split_data = split_by_region(data_file, regions)
    all_sequences = read_all_sequences(msa_file)
    results_df = build_results(split_data, regions, all_sequences, genomes, group)

    write_split_workbook(split_data, f"{group}_split.xlsx")
    results_df.to_excel(f"{group}_mutation_analysis.xlsx", index=False, sheet_name="Results")
