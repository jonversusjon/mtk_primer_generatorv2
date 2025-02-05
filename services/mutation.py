from Bio.Seq import Seq
from typing import Dict, List, Optional, Tuple
from .utils import gc_content
from .primer_select import is_overhang_compatible
from .utils import translate_codon

def find_alternative_codons(codon, codon_usage_dict, max_mutations=1, verbose=False):
    """
    Finds alternative codons for the amino acid encoded by the given codon.
    Ensures DNA-to-RNA conversion before lookup.
    """
    print(f"find_alternative_codons: {codon}, {codon_usage_dict}, {max_mutations}, {verbose}")  # Debugging
    # Convert codon to RNA format before lookup
    codon_rna = codon.replace('T', 'U')  

    # Translate codon using standard genetic code
    amino_acid = str(Seq(codon_rna).translate(table=1))

    if verbose:
        print(f"🔍 Checking codon: {codon} (RNA: {codon_rna}), Amino Acid: {amino_acid}")

    if amino_acid not in codon_usage_dict:
        print(f"⚠️ Amino acid {amino_acid} not found in codon usage dictionary.")
        return []

    # Find alternative codons
    alternative_codons = [
        {"codon": alt_codon.replace('U', 'T'), "frequency": usage}  # Convert back to DNA
        for alt_codon, usage in codon_usage_dict[amino_acid].items()
        if alt_codon != codon_rna  # Ensure original codon is excluded
    ]

    sorted_codons = sorted(alternative_codons, key=lambda x: x["frequency"], reverse=True)

    if verbose:
        print(f"✅ Alternative codons for {amino_acid}: {sorted_codons}")

    return sorted_codons


def find_codon_replacements_in_range(seq, site_start, site_end, codon_usage_dict, max_mutations, verbose=False):
    print(f"\n🔎 Searching for codon replacements in range {site_start}-{site_end}")  

    replacements = []
    seq_str = str(seq)  

    for i in range(site_start, site_end, 3):  
        codon = seq_str[i:i + 3]
        if len(codon) != 3:
            print(f"⚠️ Skipping incomplete codon {codon} at {i}")
            continue  

        # Convert codon to RNA format before lookup
        codon_rna = codon.replace('T', 'U')  

        # Translate codon to amino acid (force uppercase for dictionary lookup)
        aa = str(Seq(codon_rna).translate(table=1)).upper()
        print(f"🔍 Looking for codon: {codon} (RNA: {codon_rna}) -> Amino Acid: {aa}")  
        print(f"🛠 Codon usage dict: {codon_usage_dict}")  # Debugging
        if aa not in codon_usage_dict:
            print(f"⚠️ No codon usage data for amino acid {aa} (codon {codon})")
            print(f"🛠 Available amino acids in dict: {list(codon_usage_dict.keys())}")  # Debugging
            continue

        # Get alternative codons
        alt_codons = find_alternative_codons(codon, codon_usage_dict, max_mutations, verbose)

        if not alt_codons:
            print(f"⚠️ No alternative codons found for {codon} ({aa}) at {i}")
            continue

        replacement_entry = {
            "position": i,
            "original_codon": codon,
            "alternative_codons": alt_codons
        }

        replacements.append(replacement_entry)
        print(f"✅ Replacements for {codon} at {i}: {replacement_entry}")  

    return replacements



def gather_mutation_options(
    seq: Seq,
    sites_to_mutate: dict,
    codon_usage_dict: dict,
    spacer: str,
    bsmbi_site: str,
    min_tm: float = 57,
    template_seq: Optional[str] = None,
    verbose: bool = False,
) -> List[Tuple[Dict, Dict]]:
    """
    Designs primers with codon optimization, constraint-based design, and backtracking.
    """

    print("Gathering mutation options...")  # Debugging Step
    print(f"codon_usage_dict: {codon_usage_dict}")  # Debugging Step
    print(f"sites_to_mutate: {sites_to_mutate}")  # Debugging Step
    
    
    sites_to_mutate_list = []
    for enzyme, sites in sites_to_mutate.items():
        for site in sites:
            if "sequence" not in site or "position" not in site:
                raise ValueError(f"Site is missing required keys: {site}")

            site_sequence = str(site["sequence"])
            site_start = int(site["position"])
            site_end = site_start + len(site_sequence)
            sites_to_mutate_list.append((site_start, site_end, site_sequence))

    print(f"Sites to mutate list: {sites_to_mutate_list}")  # Debugging Step

    if not sites_to_mutate_list:
        print("⚠️ No sites to mutate found!")  # Debugging Step
        return []

    sites_to_mutate_list.sort(key=lambda x: x[0])

    # Step 1: Gather all mutation options
    mutation_sets = []
    for site_start, site_end, site_sequence in sites_to_mutate_list:
        print(f"Processing site: {site_start}-{site_end} {site_sequence}")  # Debugging Step

        mutations = get_mutations_for_site(
            seq, site_start, site_end, site_sequence, codon_usage_dict, verbose=verbose
        )

        print(f"Mutations found: {mutations}")  # Debugging Step

        # If mutations are empty, this site has no viable changes
        if not mutations:
            print(f"⚠️ No viable mutations for site {site_start}-{site_end}")  # Debugging Step
            continue

        # Ensure codon usage frequency is a number before sorting
        for mutation in mutations:
            if "codon_usage_frequency" not in mutation:
                print(f"⚠️ Mutation missing 'codon_usage_frequency': {mutation}")  # Debugging Step
                continue
            mutation["codon_usage_frequency"] = float(mutation["codon_usage_frequency"])

        mutations.sort(key=lambda x: x["codon_usage_frequency"], reverse=True)
        mutation_sets.append(mutations)

    print(f"Final mutation sets: {mutation_sets}")  # Debugging Step

    if not mutation_sets:
        print("⚠️ All mutation sets are empty!")  # Debugging Step

    return mutation_sets



def find_best_mutation_set(mutation_options, verbose=False):
    """
    Iteratively tries mutation sets, eliminating used ones to avoid retries.
    
    Args:
        mutation_options (List[List[Dict]]): Ordered list of mutation sets (mutates in place).
        verbose (bool): Whether to print debug information.
    
    Returns:
        List[Dict]: The best valid mutation set found, or an empty list if none work.
    """
    while mutation_options:
        # Get and remove the highest-priority mutation set (first in list)
        mutation_set = mutation_options.pop(0)  

        if verbose:
            print(f"\n🔍 Attempting mutation set (Remaining: {len(mutation_options)})...")

        valid_mutation_subset = find_compatible_mutation_subsets(mutation_set, verbose=verbose)

        if valid_mutation_subset:
            if verbose:
                print("✅ Found a valid mutation set!")
            return valid_mutation_subset  # Return as soon as a valid set is found

    return []  # No valid mutation set found


def get_sticky_end_options(seq, mutation_index, sticky_end_length=4):
    """
    Determines potential sticky end sequences based on a mutation location.
    """
    sticky_end_options = []

    for shift in range(-sticky_end_length + 1, 1):  
        cut_site = mutation_index + shift  

        if cut_site < 0 or cut_site + sticky_end_length > len(seq):
            continue  

        sticky_end = seq[cut_site : cut_site + sticky_end_length]

        if 0.25 <= gc_content(sticky_end) <= 0.75:  
            sticky_end_options.append({"sticky_end": sticky_end, "cut_site": cut_site})

    return sticky_end_options


def get_mutations_for_site(seq, site_start, site_end, site_sequence, codon_usage_dict, max_mutations=1, verbose=False):
    """
    Identifies mutations that disrupt restriction sites while optimizing codon usage.
    """
    print(f"\n🔍 Checking site {site_start}-{site_end} for sequence: {site_sequence}")  
    print(f"get_mutations_for_site codon_usage_dict: {codon_usage_dict}")  # Debugging
    mutations = []

    # Use find_codon_replacements_in_range to get codon-level replacements
    codon_replacements = find_codon_replacements_in_range(
        seq, site_start, site_end, codon_usage_dict, max_mutations, verbose
    )

    print(f"🧬 Codon replacements found: {codon_replacements}")  

    if not codon_replacements:
        print("⚠️ No codon replacements available for this site.")  
        return []

    for codon_data in codon_replacements:
        codon_start = codon_data["position"]
        original_codon = codon_data["original_codon"]

        for alt_codon_data in codon_data["alternative_codons"]:
            new_codon = alt_codon_data["codon"]

            for i in range(3):
                if original_codon[i] != new_codon[i]:
                    nucleotide_index = codon_start + i
                    original_nucleotide = original_codon[i]
                    new_nucleotide = new_codon[i]

                    # Apply mutation
                    mutated_seq = seq[:nucleotide_index] + new_nucleotide + seq[nucleotide_index + 1:]
                    mutated_site_region = mutated_seq[site_start:site_end]

                    if site_sequence in mutated_site_region:
                        print(f"⚠️ Rejected mutation {new_codon} at {codon_start} because {site_sequence} is still present.")
                        continue

                    # Get sticky ends
                    sticky_end_candidates = get_sticky_end_options(seq, nucleotide_index)

                    print(f"🔬 Sticky end candidates for {new_nucleotide} at {nucleotide_index}: {sticky_end_candidates}")  

                    # Validate sticky ends
                    valid_sticky_ends = [
                        se for se in sticky_end_candidates
                        if 0.25 <= gc_content(se["sticky_end"]) <= 0.75
                    ]

                    if not valid_sticky_ends:
                        print(f"⚠️ No valid sticky ends for mutation {new_codon} at {nucleotide_index}.")  
                        continue

                    # Store mutation
                    mutation_entry = {
                        "index": nucleotide_index,
                        "original_nucleotide": original_nucleotide,
                        "new_nucleotide": new_nucleotide,
                        "original_codon": original_codon,
                        "new_codon": new_codon,
                        "codon_usage_frequency": alt_codon_data["frequency"],
                        "sticky_end_options": valid_sticky_ends  
                    }
                    mutations.append(mutation_entry)
                    print(f"✅ Accepted mutation: {mutation_entry}")  

    return sorted(mutations, key=lambda x: x["codon_usage_frequency"], reverse=True)


def find_compatible_mutation_subsets(mutation_sets, verbose=False):
    """
    Finds a subset of mutations where at least one mutation per restriction site has a compatible set of sticky ends.
    """
    if verbose:
        print("\n🔬 Searching for a valid subset of point mutations...")

    all_mutation_options = []
    
    for mutation_set in mutation_sets:
        # Ensure mutation_set is always treated as a list
        if isinstance(mutation_set, dict):  # If it's a dict, wrap it in a list
            mutation_set = [mutation_set]

        mutation_options = []
        for mutation in mutation_set:  # Now mutation_set is correctly treated as a list
            if not isinstance(mutation, dict) or "sticky_end_options" not in mutation:
                print(f"⚠️ Unexpected data format: {mutation}")
                continue  # Skip malformed entries
            
            if not isinstance(mutation["sticky_end_options"], list):
                print(f"⚠️ sticky_end_options should be a list but got: {type(mutation['sticky_end_options'])}")
                continue  # Skip malformed entries
            
            for sticky_end in mutation["sticky_end_options"]:
                if not isinstance(sticky_end, dict):
                    print(f"⚠️ Expected dict but got {type(sticky_end)}: {sticky_end}")
                    continue  # Skip malformed entries
                
                mutation_options.append({
                    "mutation": mutation,
                    "sticky_end": str(sticky_end.get("sticky_end", "UNKNOWN")),  # Ensure string type
                    "cut_site": sticky_end.get("cut_site", "UNKNOWN")
                })

        all_mutation_options.append(mutation_options)

    # Try to find a valid subset
    selected_mutations = []
    for mutation_group in all_mutation_options:
        for candidate in mutation_group:
            # Check if sticky ends are compatible with previously selected ones
            if not selected_mutations or all(
                is_overhang_compatible([prev["sticky_end"], candidate["sticky_end"]])
                for prev in selected_mutations
            ):
                selected_mutations.append(candidate)
                break  # Select only one mutation per site

    if verbose:
        print(f"✅ Found valid subset with {len(selected_mutations)} compatible mutations.")

    return selected_mutations if len(selected_mutations) == len(mutation_sets) else None
