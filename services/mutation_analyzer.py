from typing import Dict, List, Optional, Tuple, Union
from Bio.Seq import Seq
from .utils import GoldenGateUtils
from .base import GoldenGateDesigner
from .mutation_optimizer import MutationOptimizer
from itertools import product
from models.mutation_dict import MutationDict

class MutationAnalyzer(GoldenGateDesigner):
    """
    Analyzes DNA sequences and generates possible mutations for Golden Gate assembly.
    Focuses on finding and evaluating potential mutation sites.
    """
    def __init__(
        self, 
        sequence: Union[Seq, str],
        codon_usage_dict: Dict[str, Dict[str, float]],
        max_mutations: int = 1,
        verbose: bool = False
    ):
        super().__init__(verbose=verbose)
        self.utils = GoldenGateUtils()
        self.mutation_optimizer = MutationOptimizer(verbose=verbose)
        self.mutation_dict = MutationDict()
        self.state = {
            'current_codon': '',
            'current_position': 0,
            'mutations_found': []
        }
        
        self.sequence = str(sequence) if isinstance(sequence, Seq) else sequence
        self.codon_usage_dict = codon_usage_dict
        self.max_mutations = max_mutations
        self.verbose = verbose

    def get_all_mutations(
        self,
        sites_to_mutate: List[Dict],
        ) -> List[List[Dict]]:
        """
        Analyze and gather all possible mutation options for given restriction sites.
        
        Returns:
            List of lists, where each inner list contains mutation options for a restriction site.
            Returns empty lists if no valid mutations are found.
        """
        try:
            all_mutations = []
            with self.debug_context("mutation_analysis"):
                for site_info in self._extract_sites_info(sites_to_mutate):
                    site_start, site_end, site_sequence, enzyme, strand = site_info
                    
                    with self.debug_context(f"site_{site_start}_{site_end}"):
                        # Step 1: Find codons that overlap with restriction sites
                        overlapping_codons = self.find_overlapping_codons(site_start, site_end)
                        self.logger.debug(f"Found {len(overlapping_codons)} overlapping codons")
                    
                    with self.debug_context(f"site_{site_start}_{site_end}_alternatives"):
                        # Step 2: Find alternative codons that disrupt the sites
                        alternative_codons = self.find_alternative_codons_for_site(overlapping_codons, site_sequence)
                        self.logger.debug(f"Found {len(alternative_codons)} alternative codons")
                    
                    with self.debug_context(f"site_{site_start}_{site_end}_mutations"):
                        # Step 3: Find disruptive mutations for the site
                        site_mutations = self._find_disruptive_mutations(site_start, site_end, site_sequence, alternative_codons, enzyme, strand)
                        self.logger.debug(f"Found {len(site_mutations)} disruptive mutations")
                    
                    if site_mutations:
                        all_mutations.append(site_mutations)

            if self.verbose:
                print(f"Found {len(all_mutations)} sites with valid mutations")
                
            return all_mutations

        except Exception as e:
            if self.verbose:
                print(f"Error in mutation analysis: {str(e)}")
            return []

    def find_overlapping_codons(self, site_start: int, site_end: int) -> List[Dict]:
        """
        Find codons that overlap with a given restriction site.
        
        Args:
            site_start: Start position of the restriction site
            site_end: End position of the restriction site
            
        Returns:
            List of dictionaries containing information about overlapping codons
        """
        self.logger.debug(f"Searching for codons overlapping with range {site_start}-{site_end}")
        overlapping_codons = []
        
        # Convert sequence to string if it's a Seq object
        sequence = str(self.sequence)
        
        # Find first codon start that could overlap with recognition site
        first_codon_start = site_start - (site_start % 3)
        if first_codon_start > site_start - 2:  # Ensure we catch partial overlap
            first_codon_start -= 3
            
        # Find last codon that could overlap
        last_codon_start = site_end - (site_end % 3)
        if last_codon_start < site_end:
            last_codon_start += 3

        for i in range(first_codon_start, last_codon_start, 3):
            codon = sequence[i:i + 3]
            if len(codon) != 3:
                continue
                
            # Determine which positions in this codon overlap with recognition site
            overlap_start = max(0, site_start - i)
            overlap_end = min(3, site_end - i)
            
            if overlap_start >= overlap_end:
                continue  # No overlap with recognition site
                
            overlapping_bases = codon[overlap_start:overlap_end]
            
            # Convert to RNA and get amino acid
            codon_rna = codon.replace('T', 'U')
            seq_obj = Seq(codon_rna)
            amino_acid = str(seq_obj.translate(table=1)).replace('*', '')
            
            overlapping_codons.append({
                "position": i,
                "codon": codon,
                "overlapping_bases": overlapping_bases,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                "aa": amino_acid
            })

        return overlapping_codons

    def find_alternative_codons_for_site(
        self, 
        overlapping_codons: List[Dict], 
        site_sequence: str,
    ) -> Dict[int, List[Dict]]:
        """
        Find alternative codons that disrupt the restriction site.
        
        Args:
            overlapping_codons: List of dictionaries containing overlapping codon information
            site_sequence: The restriction site sequence
            codon_usage_dict: Dictionary of codon usage frequencies
            max_mutations: Maximum number of mutations allowed per codon
            
        Returns:
            Dictionary mapping codon positions to lists of alternative codons
        """
        alternative_codons = {}
        
        for codon_data in overlapping_codons:
            position = codon_data["position"]
            original_codon = codon_data["codon"]
            
            # Find alternative codons
            alt_codons = self._find_alternative_codons(codon=original_codon)
            
            # Filter for codons that disrupt the site
            disruptive_alt_codons = []
            for alt in alt_codons:
                if not self._would_create_site(
                    position=position,
                    original_codon=original_codon,
                    new_codon=alt["codon"],
                    site_sequence=site_sequence,
                    overlap_start=codon_data["overlap_start"],
                    overlap_end=codon_data["overlap_end"]
                ):
                    disruptive_alt_codons.append(alt)
            
            if disruptive_alt_codons:
                alternative_codons[position] = disruptive_alt_codons
        
        return alternative_codons

    def _find_alternative_codons(
        self,
        codon: str,
    ) -> List[Dict]:
        """
        Find alternative codons for a given codon.
        
        Args:
            codon: Original codon sequence
            codon_usage_dict: Dictionary of codon usage frequencies
            max_mutations: Maximum number of mutations allowed
            
        Returns:
            List of dictionaries containing alternative codons and their frequencies
        """
        try:
            # Convert codon to RNA format
            codon_rna = codon.replace('T', 'U')
            
            # Get amino acid
            seq_obj = Seq(codon_rna)
            amino_acid = str(seq_obj.translate(table=1)).replace('*', '')
            
            self.logger.debug(f"Codon: {codon} -> RNA: {codon_rna} -> AA: {amino_acid}")
            
            # Validate dictionary
            if not isinstance(self.codon_usage_dict, dict):
                self.logger.error(f"Invalid codon usage dictionary type: {type(self.codon_usage_dict)}")
                return []
                
            if amino_acid not in self.codon_usage_dict:
                self.logger.debug(f"Amino acid {amino_acid} not found in codon usage dictionary")
                return []

            # Find alternative codons
            alternative_codons = []
            for alt_codon, usage in self.codon_usage_dict[amino_acid].items():
                if alt_codon != codon_rna:
                    dna_codon = alt_codon.replace('U', 'T')
                    num_mutations = sum(1 for i in range(3) if codon[i] != dna_codon[i])
                    if num_mutations <= self.max_mutations:
                        alternative_codons.append({
                            "codon": dna_codon,
                            "frequency": usage
                        })

            return sorted(alternative_codons, key=lambda x: x["frequency"], reverse=True)

        except Exception as e:
            self.logger.error(f"Error finding alternative codons for {codon}: {str(e)}")
            self.logger.error(f"Current state - codon: {codon}, amino acid: {amino_acid}")
            raise

    def _would_create_site(
        self,
        position: int,
        original_codon: str,
        new_codon: str,
        site_sequence: str,
        overlap_start: int,
        overlap_end: int
    ) -> bool:
        """
        Check if replacing a codon would create or maintain a restriction site.
        Uses the test sequence approach from the original working code.
        """
        sequence = str(self.sequence)
        
        # Build sequence with alternative codon
        test_seq = (
            sequence[:position] +  # Sequence before codon
            new_codon +           # Alternative codon
            sequence[position+3:] # Sequence after codon
        )
        
        # Check if site sequence is present in relevant region
        check_start = max(0, position - len(site_sequence) + 1)
        check_end = min(len(test_seq), position + len(site_sequence) + 2)
        
        region_to_check = test_seq[check_start:check_end]
        
        return site_sequence in region_to_check

    def _is_site_present(self, site_sequence: str, original_codon: str, alternative_codon: str, overlap_start: int, overlap_end: int) -> bool:
        """
        Check if the restriction site is still present after replacing the original codon with the alternative codon.
        """
        modified_site_sequence = site_sequence[:overlap_start] + alternative_codon[overlap_start:overlap_end] + site_sequence[overlap_end:]
        return site_sequence == modified_site_sequence

    def _find_disruptive_mutations(self, site_start: int, site_end: int, site_sequence: str, alternative_codons: Dict[str, List[Dict]], enzyme: str, strand: str) -> List[Dict]:
        """
        Find mutations that disrupt the restriction site.
        
        Returns:
            A list of mutation entries.
        """
        mutations = []
        
        for position, alt_codons in alternative_codons.items():
            original_codon = self.sequence[position:position+3]
            
            for alt_codon_data in alt_codons:
                alt_codon = alt_codon_data["codon"]
                frequency = alt_codon_data["frequency"]
                
                for i in range(3):
                    if original_codon[i] != alt_codon[i]:
                        mutation_entry = self._create_mutation_entry(
                            site_sequence, site_start, position, i,
                            original_codon, alt_codon, frequency,
                            enzyme, strand
                        )
                        if mutation_entry:
                            mutations.append(mutation_entry)
        
        return mutations

    def _create_mutation_entry(self, site_sequence: str, site_start: int, codon_start: int, codon_pos: int,
                               original_codon: str, new_codon: str, frequency: float, enzyme: str, strand: str) -> Optional[Dict]:
        """
        Create a mutation entry and add it to the mutation dictionary.
        
        Returns:
            The created mutation entry if successful, None otherwise.
        """
        mutation_position = codon_start + codon_pos
        original_base = original_codon[codon_pos]
        new_base = new_codon[codon_pos]
        
        # Check if the mutation disrupts the restriction site
        mutated_site_sequence = site_sequence[:codon_pos] + new_base + site_sequence[codon_pos+1:]
        if mutated_site_sequence == site_sequence:
            return None
        
        mutation_entry = {
            "position": mutation_position,
            "original_base": original_base,
            "new_base": new_base,
            "original_codon": original_codon,
            "new_codon": new_codon,
            "frequency": frequency,
            "site_start": site_start,
            "enzyme": enzyme,
            "strand": strand
        }
        
        mutation_id = f"{mutation_position}_{original_base}_{new_base}"
        self.mutation_dict.add_entry(mutation_id, mutation_entry)
        
        return mutation_entry

    def _extract_sites_info(self, sites_to_mutate: dict) -> List[Tuple[int, int, str, str, str]]:
        """Extract site information from the sites_to_mutate dictionary."""
        sites_info = []
        for enzyme, sites in sites_to_mutate.items():
            for site in sites:
                site_start = int(site["position"])
                site_end = site_start + len(site["sequence"])
                site_sequence = str(site["sequence"])
                strand = site["strand"]
                sites_info.append((site_start, site_end, site_sequence, enzyme, strand))
        return sites_info


    def process_site_mutations(
        self,
        alternative_codons: List[Dict],
        site_start: int,
        site_end: int,
        enzyme: str,
        strand: str
    ) -> List[Dict]:
        """Create and process individual mutations for a site."""
        # Implement this function to create and process mutations
        # Use the existing _create_mutation function
        # Use self.sequence instead of passing sequence as an argument
        # Return a list of dictionaries containing mutation information
        pass




































    # def get_mutable_position_bounds(self, mutation_list: List[Dict]) -> Tuple[int, int]:
    #     """
    #     Find the 5'-most and 3'-most mutable positions in a restriction site.
    #     """
    #     mutation_positions = {mut['rs_index'] for mut in mutation_list}

    #     if not mutation_positions:
    #         raise ValueError("No mutation positions found in mutation list")

    #     return min(mutation_positions), max(mutation_positions)

    # def get_valid_cut_site_range(self, five_prime_pos: int, three_prime_pos: int) -> Tuple[int, int]:
    #     """
    #     Calculate the valid range for cut sites based on mutable position bounds.
    #     """
    #     cut_start = five_prime_pos - 24  # 24bp upstream of 5'-most
    #     cut_end = three_prime_pos + 1   #  1bp downstream of 3'-most, inclusive
    #     return cut_start, cut_end
    
    # def get_sticky_end_options(self, seq: str, five_prime_pos:int, three_prime_pos:int) -> List[Dict]:
    #     """
    #     Get all valid sticky end options for a restriction site.
    #     """
    #     # five_prime_pos, three_prime_pos = self.get_mutable_position_bounds(mutation_list) # Get from mutation list
    #     cut_start, cut_end = self.get_valid_cut_site_range(five_prime_pos, three_prime_pos)
    #     self.logger.debug(f"Mutable positions: {five_prime_pos}-{three_prime_pos}")
    #     self.logger.debug(f"Valid cut site range: {cut_start}-{cut_end}")
    #     sticky_end_options = []
    #     for cut_pos in range(cut_start, cut_end):
    #         if cut_pos < 0 or cut_pos + 4 > len(seq):
    #             continue
    #         sticky_end = seq[cut_pos: cut_pos + 4]
    #         if 0.25 <= self.utils.gc_content(sticky_end) <= 0.75:
    #             sticky_end_options.append({
    #                 "sticky_end": sticky_end,
    #                 "cut_site": cut_pos
    #             })
    #     self.logger.debug(f"Found {len(sticky_end_options)} valid sticky end options")
    #     return sticky_end_options


    # def _create_mutation(
    #     self,
    #     seq: str,
    #     site_sequence: str,
    #     site_start: int,
    #     codon_start: int,
    #     position: int,
    #     original_codon: str,
    #     new_codon: str,
    #     frequency: float,
    # ) -> Optional[Dict]:
    #     """Create a single mutation entry (WITHOUT sticky end options) with debugging readouts."""

    #     nucleotide_index = codon_start + position
    #     print(f"Nucleotide Index (absolute in seq): {nucleotide_index}")

    #     # Ensure indices are within valid range
    #     if not (0 <= nucleotide_index < len(seq)):
    #         print(f"ERROR: Nucleotide index {nucleotide_index} out of bounds for sequence of length {len(seq)}")
    #         return None

    #     original_nucleotide = original_codon[position]
    #     new_nucleotide = new_codon[position]
    #     print(f"Original Nucleotide: {original_nucleotide} -> New Nucleotide: {new_nucleotide}")

    #     # Apply mutation *only to the site region*
    #     mutated_seq = seq[site_start:site_start+len(site_sequence)]
    #     mutated_nucleotide_index = nucleotide_index - site_start  # Relative index

    #     # Ensure relative index is within site sequence range
    #     if not (0 <= mutated_nucleotide_index < len(mutated_seq)):
    #         print(f"ERROR: Mutated nucleotide index {mutated_nucleotide_index} out of bounds for site sequence length {len(mutated_seq)}")
    #         return None

    #     print(f"Original Site Region: {mutated_seq}")
    #     mutated_seq = (
    #         mutated_seq[:mutated_nucleotide_index]
    #         + new_nucleotide
    #         + mutated_seq[mutated_nucleotide_index + 1:]
    #     )
    #     print(f"Mutated Site Region: {mutated_seq}")

    #     # Check if mutation actually disrupts the site
    #     if site_sequence in mutated_seq:
    #         print("Mutation failed: Restriction site still present.")
    #         return None

    #     print("Mutation successful: Restriction site disrupted.")
        
    #     mutation_entry = {
    #         "rs_index": nucleotide_index - site_start,
    #         "original_nt": original_nucleotide,
    #         "new_nt": new_nucleotide,
    #         "original_codon": original_codon,
    #         "new_codon": new_codon,
    #         "usage": frequency,
    #         "mutated_seq": mutated_seq,  # Now the 6nt (or site length) sequence
    #     }

    #     return mutation_entry


    # def _process_codon_mutations(
    #     self,
    #     seq: str,
    #     site_sequence: str,
    #     site_start: int,
    #     site_end: int,
    #     codon_data: Dict,
    #     max_mutations: int = 1,
    # ) -> List[Dict]:
    #     """Process mutations for a single codon."""
    #     mutations = []
    #     codon_start = codon_data["position"]
    #     original_codon = codon_data["original_codon"]

    #     for alt_codon_data in codon_data["alternative_codons"]:
    #         new_codon = alt_codon_data["codon"]

    #         for i in range(3):
    #             if original_codon[i] != new_codon[i]:
    #                 mutation = self._create_mutation(
    #                     seq, site_sequence, site_start,
    #                     codon_start, i, original_codon, new_codon,
    #                     alt_codon_data["frequency"]
    #                 )
    #                 if mutation:
    #                   mutation['original_codon_start'] = codon_start #add this
    #                   mutations.append(mutation)

    #     return mutations

    # def get_mutations_for_site(
    #     self,
    #     seq: str,
    #     site_start: int,
    #     site_end: int,
    #     site_sequence: str,
    #     codon_usage_dict: Dict,
    #     max_mutations: int = 1,
    # ) -> List[Dict]:
    #     """
    #     Identifies mutations that disrupt restriction sites, returns a list of
    #     dictionaries, each representing a potential mutation.  Does *not* structure
    #     the data hierarchically.  That's done in `gather_mutation_options`.
    #     """
    #     self.logger.debug(f"Analyzing site {site_start}-{site_end} for sequence: {site_sequence}")
    #     print(f"Debug - sequence type: {type(seq)}")
    #     print(f"Debug - sequence: {seq[:50]}...")  # First 50 chars
    #     print(f"Debug - site_start: {site_start}, site_end: {site_end}")
    #     print(f"Debug - site_sequence: {site_sequence}")
    #     try:
    #         seq = str(seq)
    #         site_sequence = str(site_sequence)

    #         codon_replacements = self.find_codon_replacements_in_range(
    #             seq=seq,
    #             site_start=site_start, 
    #             site_end=site_end,
    #             recognition_seq=site_sequence,
    #             codon_usage_dict=codon_usage_dict,
    #             max_mutations=max_mutations
    #         )

    #         self.logger.debug(f"Found {len(codon_replacements)} codon replacement options")

    #         if not codon_replacements:
    #             self.logger.debug("No codon replacements available for this site.")
    #             return []

    #         all_potential_mutations = []  # List of individual mutation dictionaries
    #         for codon_data in codon_replacements:
    #             new_mutations = self._process_codon_mutations(
    #                 seq, site_sequence, site_start, site_end, codon_data, max_mutations=max_mutations
    #             )
    #             all_potential_mutations.extend(new_mutations)

    #         if not all_potential_mutations:
    #             self.logger.debug("No valid mutations found after codon replacement.")
    #             return []
          
    #         return all_potential_mutations # returns a list of mutations

    #     except Exception as e:
    #         self.logger.error(f"Error getting mutations for site: {str(e)}", exc_info=True)
    #         raise
    

    
    # def _process_sites_to_mutate(self, sites_to_mutate: dict) -> List[Tuple[int, int, str, str, str]]:
    #     """Process and validate sites to mutate, including enzyme and strand."""
    #     sites_list = []
    #     for enzyme, sites in sites_to_mutate.items():
    #         for site in sites:
    #             if "sequence" not in site or "position" not in site or "strand" not in site:
    #                 raise ValueError(f"Site is missing required keys: {site}")

    #             sites_list.append((
    #                 int(site["position"]),
    #                 int(site["position"]) + len(site["sequence"]),
    #                 str(site["sequence"]),
    #                 enzyme,  # Include enzyme
    #                 site["strand"]  # Include strand: "pos" or "neg"
    #             ))
    #     return sites_list
    
    # def _get_all_mutations(
    #     self,
    #     sequence: Union[Seq, str],
    #     sites_to_mutate: List[Dict],
    #     codon_usage_dict: Dict[str, Dict[str, float]],
    #     max_mutations: int = 1,
    #     verbose: bool = False,
    # ) -> List[List[Dict]]:
    #     """
    #     Analyze and gather all possible mutation options for given restriction sites.
        
    #     Args:
    #         sequence: Input DNA sequence
    #         sites_to_mutate: List of restriction sites to mutate
    #         codon_usage_dict: Dictionary of codon usage frequencies
    #         max_mutations: Maximum number of mutations to consider
    #         verbose: Whether to print debug information
            
    #     Returns:
    #         List of lists, where each inner list contains mutation options for a restriction site.
    #         Returns empty lists if no valid mutations are found.
    #     """
    #     try:
    #         sequence_str = str(sequence) if isinstance(sequence, Seq) else sequence
    #         all_mutations = []

    #         for site_info in self._process_sites_to_mutate(sites_to_mutate):
    #             site_start, site_end, site_sequence, enzyme, strand = site_info
                
    #             # Get mutations for this specific site
    #             site_mutations = self.get_mutations_for_site(
    #                 sequence_str, site_start, site_end, site_sequence, 
    #                 codon_usage_dict, max_mutations
    #             )
                
    #             if not site_mutations:
    #                 continue
                    
    #             # Process mutations for this site
    #             site_muts = []
    #             codons_processed = set()
                
    #             for mut in site_mutations:
    #                 codon_index = mut['original_codon_start']
                    
    #                 # Create mutation entry with site information
    #                 mutation_entry = {
    #                     'rs_index': mut['rs_index'],
    #                     'original_nt': mut['original_nt'],
    #                     'new_nt': mut['new_nt'],
    #                     'new_codon': mut['new_codon'],
    #                     'usage': mut['usage'],
    #                     'mutated_seq': mut['mutated_seq'],
    #                     'site_start': site_start,
    #                     'enzyme': enzyme,
    #                     'strand': strand
    #                 }
                    
    #                 site_muts.append(mutation_entry)
                    
    #                 # Track processed codons (if needed for future reference)
    #                 codons_processed.add(codon_index)
                
    #             if site_muts:
    #                 all_mutations.append(site_muts)

    #         if verbose:
    #             print(f"Found {len(all_mutations)} sites with valid mutations")
                
    #         return all_mutations

    #     except Exception as e:
    #         if verbose:
    #             print(f"Error in mutation analysis: {str(e)}")
    #         return []
        
        
