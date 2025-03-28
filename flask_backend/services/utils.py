# services/utils.py
import csv
import json
import os
from collections import defaultdict
from functools import lru_cache
from typing import Dict, List, Optional

import numpy as np
from Bio.Data import CodonTable
from Bio.Seq import Seq
from prettytable import PrettyTable

from flask_backend.models import RestrictionSite
from flask_backend.logging import logger

class GoldenGateUtils():
    def __init__(self, verbose: bool = False, job_id: Optional[str] = None):
        self.verbose = verbose
        self.data_dir = os.path.join(
            os.path.dirname(__file__), "../static/data")
        self.codon_tables_dir = os.path.join(
            self.data_dir, "codon_usage_tables")

    def load_json_file(self, filename: str) -> Optional[Dict]:
        """Loads a JSON file from the static/data directory."""
        filepath = os.path.join(self.data_dir, filename)

        try:
            with open(filepath, "r") as file:
                return json.load(file)
        except FileNotFoundError:
            logger.error(f"File not found: {filepath}")
            return None
        except json.JSONDecodeError:
            logger.error(f"Invalid JSON format in: {filepath}")
            return None

    @lru_cache(maxsize=10)
    def get_codon_usage_dict(self, species: str) -> Optional[Dict]:
        """Loads a codon usage table for a species."""
        filename = os.path.join(self.codon_tables_dir, f"{species}.json")
        try:
            with open(filename, "r") as f:
                return json.load(f)
        except FileNotFoundError:
            logger.error(
                f"Codon usage table not found for species: {species}")
            return None

    @lru_cache(maxsize=1)
    def get_mtk_partend_sequences(self) -> Optional[Dict]:
        """Loads MTK part-end sequences."""
        return self.load_json_file("mtk_partend_sequences.json")

    def get_mtk_partend_sequence(self, mtk_part_num: str, primer_direction: str, kozak: str = "MTK") -> Optional[str]:
        """
        Retrieve the correct overhang sequence based on the part number, direction, and kozak preference.

        Args:
            mtk_part_num (str): The part number (e.g., "2", "3", "4a").
            primer_direction (str): Either "forward" or "reverse".
            kozak (str): Either "MTK" (default) or "Canonical".

        Returns:
            Optional[str]: The corresponding primer sequence if found, otherwise None.
        """
        mtk_sequences = self.get_mtk_partend_sequences()
        if not mtk_sequences:
            return None

        # Construct the standard key
        base_key = f"{mtk_part_num}{primer_direction}"

        # Check if a star version should be used for Canonical Kozak cases
        if kozak == "Canonical" and mtk_part_num in {"2", "3", "3a"} and primer_direction == "forward":
            star_key = f"{mtk_part_num}star{primer_direction}"
            return mtk_sequences.get(star_key, mtk_sequences.get(base_key))

        # Default behavior (MTK)
        return mtk_sequences.get(base_key, None)

    def get_available_species(self) -> List[str]:
        """Gets list of available species from codon usage tables."""
        if os.path.exists(self.codon_tables_dir):
            return [f[:-5] for f in os.listdir(self.codon_tables_dir)
                    if f.endswith('.json')]
        logger.warning("Codon usage tables directory not found")
        return []

    def reverse_complement(self, seq: str) -> str:
        """Returns the reverse complement of a DNA sequence."""
        return str(Seq(seq).reverse_complement())

    def get_amino_acid(self, codon: str) -> str:
        """Translates a codon to its amino acid."""
        return str(Seq(codon).translate())

    def get_codon_seqs_for_amino_acid(self, amino_acid: str) -> List[str]:
        """Returns a list of codons that encode the given amino acid."""
        print(f"Getting codons for amino acid: {amino_acid}")

        amino_acid = amino_acid.upper()
        table = CodonTable.unambiguous_dna_by_id[1]

        if amino_acid == '*':
            return table.stop_codons

        codons = [codon for codon, aa in table.forward_table.items()
                  if aa == amino_acid]
        
        print(f"Codons found: {codons}")

        return codons

    def gc_content(self, seq: str) -> float:
        """Computes GC content of a DNA sequence, rounded to 3 decimal places."""
        if not seq:
            return 0.0
        gc_count = sum(1 for nt in seq.upper() if nt in "GC")
        return round(gc_count / len(seq), 3)


    def seq_to_index(self, seq: str) -> int:
        """Converts a 4-nucleotide sequence to its corresponding matrix index."""
        seq = seq.upper()
        NT_VALUES = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        index = 0
        for pos, nt in enumerate(seq):
            power = 3 - pos
            index += NT_VALUES[nt] * (4 ** power)

        return index

    def load_compatibility_table(self, path: str) -> np.ndarray:
        """Loads the binary compatibility table into a numpy array."""
        with open(path, 'rb') as f:
            binary_data = f.read()

        compatibility_bits = np.unpackbits(
            np.frombuffer(binary_data, dtype=np.uint8)
        )
        compatibility_matrix = compatibility_bits.reshape(256, 256)

        if self.verbose:
            logger.log_step("",
                f"Loaded compatibility matrix with shape: {compatibility_matrix.shape}")
            logger.log_step("",
                f"Number of compatible pairs: {np.sum(compatibility_matrix)}")
            logger.log_step("", f"Matrix density: {np.mean(compatibility_matrix):.2%}")

        return compatibility_matrix

    def get_codon_usage(self, codon: str, amino_acid: str, codon_usage_dict: dict, default_usage: float = 0.0) -> float:
        """Retrieves codon usage frequency."""
        # Convert DNA (T) to RNA (U) for lookup
        codon_rna = codon.replace("T", "U")

        # Retrieve codon usage or default
        if amino_acid in codon_usage_dict:
            return codon_usage_dict[amino_acid].get(codon_rna, default_usage)
        return default_usage

    def convert_non_serializable(self, obj):
        """Convert non-serializable objects (Seq, ndarray, tuple) to JSON-compatible types."""
        if isinstance(obj, Seq):
            return str(obj)  # Convert BioPython Seq to string
        elif isinstance(obj, np.ndarray):
            return obj.tolist()  # Convert NumPy arrays to lists
        elif isinstance(obj, dict):
            return {k: self.convert_non_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self.convert_non_serializable(v) for v in obj]
        elif isinstance(obj, tuple):
            # Convert tuple to list for JSON
            return [self.convert_non_serializable(v) for v in obj]
        return obj  # Return as-is if it's already JSON serializable

    def analyze_incompatibility_reason(self, combo):
        """
        Analyzes why a given combination of overhangs is incompatible.
        """
        for seq in combo:
            for i in range(len(seq) - 3):
                if seq[i] == seq[i + 1] == seq[i + 2] == seq[i + 3]:
                    return f"Overhang {seq} has more than 3 consecutive identical bases."
        seen_triplets = set()
        for seq in combo:
            for i in range(len(seq) - 2):
                triplet = seq[i:i + 3]
                if triplet in seen_triplets:
                    return f"Multiple overhangs share the triplet '{triplet}'."
                seen_triplets.add(triplet)
        for seq in combo:
            gc_content = sum(
                1 for base in seq if base in "GC") / len(seq) * 100
            if gc_content == 0:
                return f"Overhang {seq} has 0% GC content (all A/T)."
            elif gc_content == 100:
                return f"Overhang {seq} has 100% GC content (all G/C)."
        return "Unknown reason (should not happen)."

    def calculate_optimal_primer_length(self, sequence: str, position: int, direction='forward') -> int:
        logger.log_step("Calculate Primer Length", f"Optimal {direction} primer from pos {position}",
                      {"sequence_length": len(sequence)})
        min_length, max_length, target_tm = 18, 30, 60
        optimal_length = min_length

        if direction == 'forward':
            for length in range(min_length, min(max_length + 1, len(sequence) - position)):
                primer_seq = sequence[position:position + length]
                tm = self.calculate_tm(primer_seq)
                logger.log_step("Length Iteration", f"Length {length}", {
                              "tm": tm, "target": target_tm})
                if tm >= target_tm:
                    optimal_length = length
                    break
        else:
            for length in range(min_length, min(max_length + 1, position + 1)):
                if position - length < 0:
                    break
                primer_seq = sequence[position - length:position]
                tm = self.calculate_tm(primer_seq)
                logger.log_step("Length Iteration", f"Length {length}", {
                              "tm": tm, "target": target_tm})
                if tm >= target_tm:
                    optimal_length = length
                    break

        logger.validate(
            optimal_length >= min_length,
            f"Calculated optimal primer length: {optimal_length}",
            {"direction": direction, "min_length": min_length, "max_length": max_length}
        )
        return optimal_length

    def calculate_tm(self, sequence: str) -> float:
        if not sequence:
            return 0.0
        sequence = sequence.upper()
        length = len(sequence)
        a_count = sequence.count("A")
        t_count = sequence.count("T")
        g_count = sequence.count("G")
        c_count = sequence.count("C")
        if length < 14:
            tm = (a_count + t_count) * 2 + (g_count + c_count) * 4
        else:
            tm = 64.9 + (41 * (g_count + c_count - 16.4)) / length
        return round(tm, 2)

    def export_primers_to_tsv(
        self,
        output_tsv_path: str,
        primer_data: Optional[List[List[str]]] = None,
        forward_primers: Optional[List[tuple]] = None,
        reverse_primers: Optional[List[tuple]] = None,
        header: Optional[List[str]] = None
    ) -> None:
        """
        Exports primers to a TSV file.

        The function works in one of two modes:
        1. If `primer_data` is provided, it writes that data as rows (assuming each row is a list of strings).
        2. If `primer_data` is None but both `forward_primers` and `reverse_primers` are provided,
            it writes both lists as rows with a default message in the third column.

        Optionally, a custom header can be provided; otherwise, a default header is used.
        """
        # Define a default header if none is provided
        if header is None:
            header = ["Primer Name", "Sequence", "Amplicon"]

        try:
            with open(output_tsv_path, mode="w", newline="") as tsv_file:
                writer = csv.writer(tsv_file, delimiter="\t")
                writer.writerow(header)

                # If primer_data is provided, write it directly
                if primer_data is not None:
                    if not primer_data:
                        logger.warning("No primer data to save.")
                        return
                    for row in primer_data:
                        writer.writerow(list(map(str, row)))
                # Otherwise, if forward and reverse primers are provided, merge them.
                elif forward_primers is not None and reverse_primers is not None:
                    # Define the default message for assembly
                    assembly_message = "Generated for Golden Gate Assembly"
                    for name, sequence in forward_primers:
                        writer.writerow([name, sequence, assembly_message])
                    for name, sequence in reverse_primers:
                        writer.writerow([name, sequence, assembly_message])
                else:
                    logger.error("Insufficient primer data provided.")
                    raise ValueError(
                        "Either primer_data or both forward_primers and reverse_primers must be provided.")

            logger.log_step("", f"Primers exported to {output_tsv_path}")
        except IOError as e:
            logger.error(f"Error writing to file {output_tsv_path}: {e}")
            raise

    def get_nested_keys(self, item, prefix='', depth=0, max_depth=100):
        """Recursively collect keys from dictionaries and lists."""
        if depth > max_depth:
            return [f"Max recursion depth ({max_depth}) exceeded."]
        keys = []
        if isinstance(item, dict):
            for k, v in item.items():
                full_key = f"{prefix}.s{k}" if prefix else k
                keys.append(full_key)
                keys.extend(self.get_nested_keys(
                    v, full_key, depth + 1, max_depth))
        elif isinstance(item, list):
            for index, element in enumerate(item):
                indexed_prefix = f"{prefix}[{index}]"
                keys.extend(self.get_nested_keys(
                    element, indexed_prefix, depth + 1, max_depth))
        return keys

    def print_object_schema(self, obj, indent=0, name="object"):
        """
        Recursively prints a text-based schema of an object.
        For dictionaries and objects (via __dict__), if multiple keys/attributes have the same
        structure, only one representative is printed along with a note of additional items.

        Parameters:
        obj: The object to schematize.
        indent: Current indentation level.
        name: Name/label for the current object.
        """
        prefix = "  " * indent

        if isinstance(obj, dict):
            print(f"{prefix}{name} (dict) with {len(obj)} keys:")
            # Group keys by the structure of their corresponding value
            groups = {}
            for key, value in obj.items():
                sig = self.get_structure(value)
                groups.setdefault(sig, []).append((key, value))
            # Print one representative per group
            for sig, items in groups.items():
                # Sort keys for consistency
                items.sort(key=lambda x: str(x[0]))
                rep_key, rep_val = items[0]
                key_repr = f"['{rep_key}']" if isinstance(
                    rep_key, str) else f"[{rep_key}]"
                self.print_object_schema(
                    rep_val, indent + 1, f"{name}{key_repr}")
                if len(items) > 1:
                    print(
                        f"{prefix}  ... ({len(items)-1} more keys with same schema)")

        elif isinstance(obj, list):
            print(f"{prefix}{name} (list) with {len(obj)} items:")
            if not obj:
                return
            # Check if all list items share the same schema
            first_sig = self.get_structure(obj[0])
            all_same = all(self.get_structure(item) ==
                           first_sig for item in obj)
            if all_same:
                self.print_object_schema(obj[0], indent + 1, f"{name}[0]")
                if len(obj) > 1:
                    print(
                        f"{prefix}  ... ({len(obj)-1} more items with same schema)")
            else:
                # Otherwise, print the schema for each item
                for idx, item in enumerate(obj):
                    self.print_object_schema(
                        item, indent + 1, f"{name}[{idx}]")

        elif hasattr(obj, "__dict__") and obj.__dict__:
            print(f"{prefix}{name} (object):")
            # Group object attributes by their schema
            groups = {}
            for attr, value in obj.__dict__.items():
                sig = self.get_structure(value)
                groups.setdefault(sig, []).append((attr, value))
            for sig, items in groups.items():
                items.sort(key=lambda x: str(x[0]))
                rep_attr, rep_val = items[0]
                self.print_object_schema(
                    rep_val, indent + 1, f"{name}.{rep_attr}")
                if len(items) > 1:
                    print(
                        f"{prefix}  ... ({len(items)-1} more attributes with same schema)")

        elif hasattr(obj, "shape") and hasattr(obj, "dtype"):
            # For array-like objects (e.g., numpy arrays, pandas DataFrames)
            print(
                f"{prefix}{name} (array-like): shape={obj.shape}, dtype={obj.dtype}")

        else:
            # Primitive or other types
            val_str = repr(obj)
            display_val = f"{val_str[:50]}{'...' if len(val_str) > 50 else ''}"
            print(f"{prefix}{name} ({type(obj).__name__}): {display_val}")

    def get_structure(self, obj):
        """
        Recursively computes a hashable signature representing the schema of `obj`.
        This signature can be used to compare whether two objects share the same structure.

        Returns a tuple describing the structure.
        """
        if isinstance(obj, dict):
            # For dictionaries, use a sorted tuple of (key, schema) pairs.
            structure = tuple(sorted(
                ((str(key), self.get_structure(value))
                 for key, value in obj.items())
            ))
            return ("dict", structure)
        elif isinstance(obj, list):
            if not obj:
                return ("list", "empty")
            # Compute the structure for each item.
            item_structs = [self.get_structure(item) for item in obj]
            # If all items share the same structure, collapse to one signature.
            if all(s == item_structs[0] for s in item_structs):
                return ("list", item_structs[0])
            else:
                return ("list", tuple(item_structs))
        elif hasattr(obj, "__dict__") and obj.__dict__:
            structure = tuple(sorted(
                ((attr, self.get_structure(value))
                 for attr, value in obj.__dict__.items())
            ))
            return ("object", structure)
        elif hasattr(obj, "shape") and hasattr(obj, "dtype"):
            return ("array", str(obj.shape), str(obj.dtype))
        else:
            return type(obj).__name__

    def summarize_bsmbi_bsai_sites(self, restriction_sites: List[RestrictionSite]) -> None:
        """
        Creates a formatted summary of restriction sites.
        """
        site_type_descriptions = {
            'BsmBI': 'BsmBI Restriction Site',
            'BsaI': 'BsaI Restriction Site',
        }

        # Group the sites by enzyme
        grouped_sites = defaultdict(list)
        for site in restriction_sites:
            grouped_sites[site.enzyme].append(site)

        table = PrettyTable()
        table.field_names = ["Site Type", "Number of Instances", "Position(s)"]

        for enzyme, sites in grouped_sites.items():
            if not sites:
                continue

            enzyme_desc = site_type_descriptions.get(enzyme, enzyme)
            positions = ", ".join(str(site.position + 1) for site in sites)
            table.add_row([enzyme_desc, len(sites), positions])

        logger.log_step("", "\nRestriction Site Analysis Summary:")
        logger.log_step("", f"\n{table}")