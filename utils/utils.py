import json
import os
import csv
from Bio.Seq import Seq
from Bio import SeqIO
from prettytable import PrettyTable

def get_codon_usage_table(species, data_dir="data/codon_usage_tables"):
    """
    Loads a codon usage table from a JSON file.

    Args:
        species (str): The name of the species (e.g., "homo_sapiens").
        data_dir (str): The directory containing the codon usage tables.

    Returns:
        dict: The codon usage table as a nested dictionary, or None if the file is not found.
    """
    filename = os.path.join(data_dir, f"{species}.json")
    try:
        with open(filename, "r") as f:
            codon_usage_table = json.load(f)
        return codon_usage_table
    except FileNotFoundError:
        return None

def get_mtk_partend_sequences(data_dir="data"):
    """
    Loads MTK part-end sequences from a JSON file.

    Args:
        data_dir (str): The directory containing the JSON file.

    Returns:
        dict: The MTK part-end sequences, or None if the file is not found.
    """
    filename = os.path.join(data_dir, "mtk_partend_sequences.json")
    try:
        with open(filename, "r") as f:
            partend_sequences = json.load(f)
        return partend_sequences
    except FileNotFoundError:
        return None


def get_amino_acid(codon):
    return Seq(codon).translate()


def get_available_species(data_dir="data/codon_usage_tables"):
    """
    Gets a list of available species from JSON files in the specified data directory.

    Parameters:
        data_dir (str): Path to the directory containing codon usage JSON files.

    Returns:
        list[str]: A list of species names derived from the JSON filenames.
                   Returns an empty list if no files are found or if the directory does not exist.
    """
    if os.path.exists(data_dir):
        return [f[:-5] for f in os.listdir(data_dir) if f.endswith('.json')]
    return []  # Return an empty list if the directory is missing


def export_primers_to_tsv(forward_primers, reverse_primers, filename="primers.tsv"):
    """
    Exports the forward and reverse primers to a TSV file.

    Args:
        forward_primers (list): List of tuples containing forward primer names and sequences.
        reverse_primers (list): List of tuples containing reverse primer names and sequences.
        filename (str): Name of the TSV file to create.
    """
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        
        for name, sequence in forward_primers:
            writer.writerow([name, sequence, "Generated for Golden Gate Assembly"])
        
        for name, sequence in reverse_primers:
            writer.writerow([name, sequence, "Generated for Golden Gate Assembly"])


def process_uploaded_file(uploaded_file):
    """Processes the uploaded plasmid file to extract sequence, annotations, and other data."""
    if not uploaded_file:
        raise ValueError("No file uploaded.")

    # Handle different possible structures
    if isinstance(uploaded_file, dict):
        files = uploaded_file.items()  # Dict structure
    elif isinstance(uploaded_file, (list, tuple)):
        files = enumerate(uploaded_file)  # Tuple or list structure
    else:
        raise ValueError(f"Unsupported uploaded file format: {type(uploaded_file)}")

    results = []
    for file in files:
        if isinstance(file, tuple) and len(file) == 2:
            filename, file_metadata = file
        elif isinstance(file, dict):
            file_metadata = file
            filename = file_metadata.get("name", "unknown_file")
        else:
            raise ValueError(f"Unrecognized file structure: {file}")

        file_content = file_metadata.get("content")
        if not file_content:
            raise ValueError(f"No content found for file: {filename}")

        _, ext = os.path.splitext(filename)

        # Process file by type
        if ext.lower() == ".dna":
            from snapgene_reader import snapgene_file_to_dict
            data = snapgene_file_to_dict(file_content)
            results.append({
                "filename": filename,
                "sequence": data["seq"],
                "annotations": data["features"],
                "primers": data.get("primers", []),
            })
        elif ext.lower() in {".gb", ".gbk"}:
            from Bio import SeqIO
            from io import BytesIO
            record = SeqIO.read(BytesIO(file_content), "genbank")
            results.append({
                "filename": filename,
                "sequence": str(record.seq),
                "annotations": record.features,
                "primers": [],
            })
        elif ext.lower() == ".ape":
            content_str = file_content.decode("utf-8")
            sequence = content_str.split("ORIGIN")[-1].replace("\n", "").replace(" ", "")
            results.append({
                "filename": filename,
                "sequence": sequence,
                "annotations": [],
                "primers": [],
            })
        else:
            raise ValueError(f"Unsupported file type: {ext}")

    return results


def plot_primers_for_mutation_site(seq, primers, output_file="dna_diagram.pdf"):
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.Graphics import GenomeDiagram
    from reportlab.lib import colors
    import os

    gd_diagram = GenomeDiagram.Diagram("DNA Diagram")
    gd_track_for_features = gd_diagram.new_track(1, name="Features")
    gd_feature_set = gd_track_for_features.new_set()

    # Add primers as features
    for primer in primers:
        primer_name = primer["primer_name"]
        start = primer["mutation_position"]
        end = start + len(primer["primer_sequence"])
        strand = 1 if "_FW" in primer_name else -1
        color = colors.green if strand == 1 else colors.red

        # Use FeatureLocation to specify strand and location
        feature_location = FeatureLocation(start, end, strand=strand)
        gd_feature_set.add_feature(
            SeqFeature(location=feature_location),
            color=color,
            label=True,
            name=primer_name
        )

    try:
        gd_diagram.draw(format="linear", pagesize="A4", fragments=1, start=0, end=len(seq))
        gd_diagram.write(output_file, "PDF")
        print(f"PDF created successfully: {os.path.abspath(output_file)}")
    except Exception as e:
        print(f"Error creating PDF: {e}")

    # Verify that the file exists after writing
    if os.path.exists(output_file):
        print(f"File saved successfully at: {os.path.abspath(output_file)}")
    else:
        print("Error: PDF file was not created.")


def gc_content(seq):
    """
    Computes the GC content of a DNA sequence as a fraction.
    """
    if len(seq) == 0:
        return 0
    return sum(1 for nt in seq if nt in "GC") / len(seq)


from itertools import product
from prettytable import PrettyTable

def rank_and_print_mutation_sets(mutation_options, verbose):
    """
    Generates all possible mutation sets, ranks them by total codon usage frequency,
    and prints them in a PrettyTable format.

    Args:
        mutation_options (List[List[Dict]]): A list of lists, where each inner list contains
                                             mutations for one restriction site.
        verbose (bool): Whether to print details.

    Returns:
        List[Tuple[List[Dict], float]]: Sorted list of mutation sets with their total codon usage score.
    """
    ranked_mutation_sets = []

    # Generate all possible mutation combinations (Cartesian product of mutation choices)
    for mutation_set in product(*mutation_options):
        total_score = sum(mutation["codon_usage_frequency"] for mutation in mutation_set)
        ranked_mutation_sets.append((mutation_set, total_score))

    # Sort by total codon usage frequency (higher is better)
    ranked_mutation_sets.sort(key=lambda x: x[1], reverse=True)

    # Print ranked mutation sets using PrettyTable
    if verbose:
        print("\n🔬 **Ranked Mutation Sets (Highest Codon Usage First)** 🔬\n")
        for rank, (mutations, score) in enumerate(ranked_mutation_sets, 1):
            table = PrettyTable()
            table.title = f"Mutation Set {rank} (Total Codon Usage Score: {score:.2f})"
            table.field_names = ["Site #", "Original Codon", "New Codon", "Usage Frequency (%)"]
            
            for site_index, mut in enumerate(mutations, start=1):
                table.add_row([site_index, mut["original_codon"], mut["new_codon"], f"{mut['codon_usage_frequency']:.1f}"])
            
            print(table)
            print("---")

    return ranked_mutation_sets  # Return sorted mutation sets
