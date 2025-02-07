from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq
from .utils import GoldenGateUtils
from .base import GoldenGateDesigner
from .mutation_optimizer import MutationOptimizer
from itertools import product

class MutationAnalyzer(GoldenGateDesigner):
    """
    Analyzes DNA sequences and generates possible mutations for Golden Gate assembly.
    Focuses on finding and evaluating potential mutation sites.
    """
    def __init__(self, verbose: bool = False):
        super().__init__(verbose=verbose)  # Call superclass constructor
        self.utils = GoldenGateUtils()
        self.mutation_optimizer = MutationOptimizer(verbose=verbose)
        self.state = {
            'current_codon': '',
            'current_position': 0,
            'mutations_found': []
        }

    def find_alternative_codons(self, codon: str, codon_usage_dict: Dict, max_mutations: int = 1) -> List[Dict]:
        """
        Finds alternative codons for the amino acid encoded by the given codon.
        """
        try:
            # Convert codon to RNA format before lookup
            codon_rna = codon.replace('T', 'U')
            
            # Create Seq object
            seq_obj = Seq(codon_rna)
            # Get amino acid and force to string, remove any asterisk for stop codons
            amino_acid = str(seq_obj.translate(table=1)).replace('*', '')
            
            self.logger.debug(f"Codon: {codon} -> RNA: {codon_rna} -> AA: {amino_acid}")
            
            # Check dictionary
            if not isinstance(codon_usage_dict, dict):
                self.logger.error(f"Invalid codon usage dictionary type: {type(codon_usage_dict)}")
                return []
                
            if amino_acid not in codon_usage_dict:
                self.logger.debug(f"Amino acid {amino_acid} not found in codon usage dictionary")
                return []

            # Find alternative codons
            alternative_codons = []
            for alt_codon, usage in codon_usage_dict[amino_acid].items():
                if alt_codon != codon_rna:
                    dna_codon = alt_codon.replace('U', 'T')
                    num_mutations = sum(1 for i in range(3) if codon[i] != dna_codon[i])
                    if num_mutations <= max_mutations:
                        alternative_codons.append({
                            "codon": dna_codon,
                            "frequency": usage
                        })

            return sorted(alternative_codons, key=lambda x: x["frequency"], reverse=True)

        except Exception as e:
            self.logger.error(f"Error finding alternative codons for {codon}: {str(e)}")
            self.logger.error(f"Current state - codon: {codon}, translation result type: {type(amino_acid)}")
            raise

    def find_codon_replacements_in_range(
        self,
        seq: str,
        site_start: int,
        site_end: int,
        recognition_seq: str,  # Add recognition sequence as parameter
        codon_usage_dict: Dict,
        max_mutations: int = 1
    ) -> List[Dict]:
        """Find codon replacements within a specified range."""
        self.logger.debug(f"Searching for codon replacements in range {site_start}-{site_end}")
        replacements = []
        seq_str = str(seq)
        
        # Find first codon start that could overlap with recognition site
        first_codon_start = site_start - (site_start % 3)
        if first_codon_start > site_start - 2:  # Ensure we catch partial overlap
            first_codon_start -= 3
            
        # Find last codon that could overlap
        last_codon_start = site_end - (site_end % 3)
        if last_codon_start < site_end:
            last_codon_start += 3

        for i in range(first_codon_start, last_codon_start, 3):
            codon = seq_str[i:i + 3]
            if len(codon) != 3:
                continue
                
            # Determine which positions in this codon overlap with recognition site
            codon_overlap_start = max(0, site_start - i)
            codon_overlap_end = min(3, site_end - i)
            
            if codon_overlap_start >= codon_overlap_end:
                continue  # No overlap with recognition site
                
            alt_codons = self.find_alternative_codons(codon, codon_usage_dict, max_mutations)
            valid_alternatives = []
            
            for alt in alt_codons:
                alt_codon = alt["codon"]
                # Build sequence with alternative codon to check if it disrupts recognition site
                test_seq = (
                    seq_str[:i] +  # Sequence before codon
                    alt_codon +    # Alternative codon
                    seq_str[i+3:]  # Sequence after codon
                )
                
                # Check if recognition sequence is disrupted
                if recognition_seq not in test_seq[site_start:site_end]:
                    valid_alternatives.append(alt)
            
            if valid_alternatives:
                replacement_entry = {
                    "position": i,
                    "original_codon": codon,
                    "alternative_codons": valid_alternatives,
                    "recognition_site_overlap": {
                        "start": codon_overlap_start,
                        "end": codon_overlap_end
                    }
                }
                replacements.append(replacement_entry)

        return replacements

    def get_mutable_position_bounds(self, mutation_list: List[Dict]) -> Tuple[int, int]:
        """
        Find the 5'-most and 3'-most mutable positions in a restriction site.
        """
        mutation_positions = {mut['rs_index'] for mut in mutation_list}

        if not mutation_positions:
            raise ValueError("No mutation positions found in mutation list")

        return min(mutation_positions), max(mutation_positions)

    def get_valid_cut_site_range(self, five_prime_pos: int, three_prime_pos: int) -> Tuple[int, int]:
        """
        Calculate the valid range for cut sites based on mutable position bounds.
        """
        cut_start = five_prime_pos - 24  # 24bp upstream of 5'-most
        cut_end = three_prime_pos + 1   #  1bp downstream of 3'-most, inclusive
        return cut_start, cut_end
    
    def get_sticky_end_options(self, seq: str, five_prime_pos:int, three_prime_pos:int) -> List[Dict]:
        """
        Get all valid sticky end options for a restriction site.
        """
        # five_prime_pos, three_prime_pos = self.get_mutable_position_bounds(mutation_list) # Get from mutation list
        cut_start, cut_end = self.get_valid_cut_site_range(five_prime_pos, three_prime_pos)
        self.logger.debug(f"Mutable positions: {five_prime_pos}-{three_prime_pos}")
        self.logger.debug(f"Valid cut site range: {cut_start}-{cut_end}")
        sticky_end_options = []
        for cut_pos in range(cut_start, cut_end):
            if cut_pos < 0 or cut_pos + 4 > len(seq):
                continue
            sticky_end = seq[cut_pos: cut_pos + 4]
            if 0.25 <= self.utils.gc_content(sticky_end) <= 0.75:
                sticky_end_options.append({
                    "sticky_end": sticky_end,
                    "cut_site": cut_pos
                })
        self.logger.debug(f"Found {len(sticky_end_options)} valid sticky end options")
        return sticky_end_options


    def _create_mutation(
        self,
        seq: str,
        site_sequence: str,
        site_start: int,
        codon_start: int,
        position: int,
        original_codon: str,
        new_codon: str,
        frequency: float,
    ) -> Optional[Dict]:
        """Create a single mutation entry (WITHOUT sticky end options) with debugging readouts."""

        nucleotide_index = codon_start + position
        print(f"Nucleotide Index (absolute in seq): {nucleotide_index}")

        # Ensure indices are within valid range
        if not (0 <= nucleotide_index < len(seq)):
            print(f"ERROR: Nucleotide index {nucleotide_index} out of bounds for sequence of length {len(seq)}")
            return None

        original_nucleotide = original_codon[position]
        new_nucleotide = new_codon[position]
        print(f"Original Nucleotide: {original_nucleotide} -> New Nucleotide: {new_nucleotide}")

        # Apply mutation *only to the site region*
        mutated_seq = seq[site_start:site_start+len(site_sequence)]
        mutated_nucleotide_index = nucleotide_index - site_start  # Relative index

        # Ensure relative index is within site sequence range
        if not (0 <= mutated_nucleotide_index < len(mutated_seq)):
            print(f"ERROR: Mutated nucleotide index {mutated_nucleotide_index} out of bounds for site sequence length {len(mutated_seq)}")
            return None

        print(f"Original Site Region: {mutated_seq}")
        mutated_seq = (
            mutated_seq[:mutated_nucleotide_index]
            + new_nucleotide
            + mutated_seq[mutated_nucleotide_index + 1:]
        )
        print(f"Mutated Site Region: {mutated_seq}")

        # Check if mutation actually disrupts the site
        if site_sequence in mutated_seq:
            print("Mutation failed: Restriction site still present.")
            return None

        print("Mutation successful: Restriction site disrupted.")
        
        mutation_entry = {
            "rs_index": nucleotide_index - site_start,
            "original_nt": original_nucleotide,
            "new_nt": new_nucleotide,
            "original_codon": original_codon,
            "new_codon": new_codon,
            "usage": frequency,
            "mutated_seq": mutated_seq,  # Now the 6nt (or site length) sequence
        }

        return mutation_entry


    def _process_codon_mutations(
        self,
        seq: str,
        site_sequence: str,
        site_start: int,
        site_end: int,
        codon_data: Dict,
        max_mutations: int = 1,
    ) -> List[Dict]:
        """Process mutations for a single codon."""
        mutations = []
        codon_start = codon_data["position"]
        original_codon = codon_data["original_codon"]

        for alt_codon_data in codon_data["alternative_codons"]:
            new_codon = alt_codon_data["codon"]

            for i in range(3):
                if original_codon[i] != new_codon[i]:
                    mutation = self._create_mutation(
                        seq, site_sequence, site_start,
                        codon_start, i, original_codon, new_codon,
                        alt_codon_data["frequency"]
                    )
                    if mutation:
                      mutation['original_codon_start'] = codon_start #add this
                      mutations.append(mutation)

        return mutations

    def get_mutations_for_site(
        self,
        seq: str,
        site_start: int,
        site_end: int,
        site_sequence: str,
        codon_usage_dict: Dict,
        max_mutations: int = 1,
    ) -> List[Dict]:
        """
        Identifies mutations that disrupt restriction sites, returns a list of
        dictionaries, each representing a potential mutation.  Does *not* structure
        the data hierarchically.  That's done in `gather_mutation_options`.
        """
        self.logger.debug(f"Analyzing site {site_start}-{site_end} for sequence: {site_sequence}")
        print(f"Debug - sequence type: {type(seq)}")
        print(f"Debug - sequence: {seq[:50]}...")  # First 50 chars
        print(f"Debug - site_start: {site_start}, site_end: {site_end}")
        print(f"Debug - site_sequence: {site_sequence}")
        try:
            seq = str(seq)
            site_sequence = str(site_sequence)

            codon_replacements = self.find_codon_replacements_in_range(
                seq=seq,
                site_start=site_start, 
                site_end=site_end,
                recognition_seq=site_sequence,
                codon_usage_dict=codon_usage_dict,
                max_mutations=max_mutations
            )

            self.logger.debug(f"Found {len(codon_replacements)} codon replacement options")

            if not codon_replacements:
                self.logger.debug("No codon replacements available for this site.")
                return []

            all_potential_mutations = []  # List of individual mutation dictionaries
            for codon_data in codon_replacements:
                new_mutations = self._process_codon_mutations(
                    seq, site_sequence, site_start, site_end, codon_data, max_mutations=max_mutations
                )
                all_potential_mutations.extend(new_mutations)

            if not all_potential_mutations:
                self.logger.debug("No valid mutations found after codon replacement.")
                return []
          
            return all_potential_mutations # returns a list of mutations

        except Exception as e:
            self.logger.error(f"Error getting mutations for site: {str(e)}", exc_info=True)
            raise
    
    def gather_mutation_options(
        self,
        seq: Seq,
        sites_to_mutate: Dict,
        codon_usage_dict: Dict,
        max_mutations: int,
        verbose: bool
    ) -> Dict:
        """Gathers all possible mutation options, organized by restriction site."""

        mutation_data = {"restriction_sites": []}

        for site_info in self._process_sites_to_mutate(sites_to_mutate):
            site_start, site_end, site_sequence, enzyme, strand = site_info
            site_mutations = self.get_mutations_for_site(
                str(seq), site_start, site_end, site_sequence, codon_usage_dict, max_mutations
            )

            if not site_mutations:  # Skip sites with no valid mutations
                continue

            # Extract codons and their mutations, restructure
            codons_data = {}
            for mut in site_mutations:
                codon_index = mut['original_codon_start']  # Use the original codon start
                if codon_index not in codons_data:
                    codons_data[codon_index] = {
                        "codon_index": codon_index,
                        "original_codon": mut["original_codon"],
                        "mutations": [],
                    }
                # Remove 'mutated_sequence' and sticky ends from *individual* mutation
                codons_data[codon_index]["mutations"].append({
                    "rs_index": mut["rs_index"],
                    "original_nt": mut["original_nt"],
                    "new_nt": mut["new_nt"],
                    "new_codon": mut["new_codon"],
                    "usage": mut["usage"],
                     "mutated_seq": mut["mutated_seq"]
                })
            # Find the best mutated_sequence to determine the best sticky_ends:
            best_mutated_sequence = ""
            best_freq = 0
            for mutation in site_mutations:
                if mutation["usage"] > best_freq:
                    best_mutated_sequence = mutation["mutated_seq"]
                    best_freq = mutation["usage"]

            #  Calculate sticky ends *once per site*
            five_prime_relative, three_prime_relative = self.get_mutable_position_bounds(site_mutations)
            # sticky_end_options = self.get_sticky_end_options(best_mutated_sequence, five_prime_relative, three_prime_relative)


            mutation_data["restriction_sites"].append({
                "enzyme": enzyme,
                "strand": strand,
                "start_index": site_start,
                "codons": list(codons_data.values()),  # Convert to list
                # "sticky_end_options": sticky_end_options
            })
            print(f"mutation_data: {mutation_data}")
        return mutation_data
    
    def _process_sites_to_mutate(self, sites_to_mutate: dict) -> List[Tuple[int, int, str, str, str]]:
        """Process and validate sites to mutate, including enzyme and strand."""
        sites_list = []
        for enzyme, sites in sites_to_mutate.items():
            for site in sites:
                if "sequence" not in site or "position" not in site or "strand" not in site:
                    raise ValueError(f"Site is missing required keys: {site}")

                sites_list.append((
                    int(site["position"]),
                    int(site["position"]) + len(site["sequence"]),
                    str(site["sequence"]),
                    enzyme,  # Include enzyme
                    site["strand"]  # Include strand: "pos" or "neg"
                ))
        return sites_list

    def find_best_mutation_set(self, mutation_options: List[Dict], verbose: bool) -> List[Dict]:
        """
        Selects the best set of mutations based on codon usage and sticky end compatibility.
        This is a simplified example and needs more sophisticated logic for real-world use.
        """

        if not mutation_options:
            return []
        
        #Simple implementation, selects mutations with the best codon usage for each index.
        best_mutations_by_index = {}
        for mutation in mutation_options:
            index = mutation["rs_index"]
            if index not in best_mutations_by_index or mutation["usage"] > best_mutations_by_index[index]["usage"]:
                best_mutations_by_index[index] = mutation

        selected_mutations = list(best_mutations_by_index.values())

        if verbose:
            self.logger.info(f"Selected mutations: {selected_mutations}")

        return selected_mutations