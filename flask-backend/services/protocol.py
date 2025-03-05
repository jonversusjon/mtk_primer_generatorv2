# services/protocol.py

from typing import List, Dict, Optional
from .sequence_prep import SequencePreparator
from .primer_design import PrimerDesigner
from .primer_select import PrimerSelector
from .mutation_analyzer import MutationAnalyzer
from .mutation_optimizer import MutationOptimizer
from .utils import GoldenGateUtils
from config.logging_config import logger
from .base import debug_context


class GoldenGateProtocol:
    """
    Orchestrates the Golden Gate protocol by managing sequence preparation,
    primer design, mutation analysis, and optimization.
    """

    def __init__(
        self,
        sequencesToDomesticate: List[Dict[str, str]],
        codon_usage_dict: Dict[str, Dict[str, float]],
        max_mutations: int,
        template_seq: Optional[str] = None,
        kozak: str = "MTK",
        output_tsv_path: str = "designed_primers.tsv",
        verbose: bool = False
    ):
        self.logger = logger.getChild("GoldenGateProtocol")
        self.utils = GoldenGateUtils()
        self.sequence_preparator = SequencePreparator()
        self.primer_designer = PrimerDesigner(kozak=kozak, verbose=verbose)
        self.primer_selector = PrimerSelector()
        self.mutation_optimizer = MutationOptimizer(verbose=verbose)
        self.logger.debug(
            f"GoldenGateProtocol initialized with codon_usage_dict: {codon_usage_dict}")
        if verbose:
            self.logger.info("GoldenGateProtocol is running in verbose mode.")
        self.mutation_analyzer = MutationAnalyzer(
            codon_usage_dict=codon_usage_dict,
            max_mutations=max_mutations,
            verbose=verbose
        )

        self.sequencesToDomesticate = sequencesToDomesticate
        self.template_seq = template_seq
        self.kozak = kozak
        self.verbose = verbose
        self.codon_usage_dict = codon_usage_dict
        self.max_mutations = max_mutations
        self.output_tsv_path = output_tsv_path

        self.mtk_partend_sequences = self.utils.get_mtk_partend_sequences()
        self.state = {
            'current_sequence_index': 0,
            'current_step': '',
            'mutations_found': [],
            'primers_designed': []
        }

    def create_gg_protocol(self) -> dict:
        """
        Main function to orchestrate the Golden Gate protocol creation.
        Returns:
            dict: A dictionary containing protocol details.
        """
        logger.info("Starting Golden Gate protocol creation...")

        result_data = {}

        for i, seq_object in enumerate(self.sequencesToDomesticate):
            sequence_number = str(i + 1)
            single_seq = seq_object["sequence"]

            sequence_data = {
                "sequence_index": i,
                "processed_sequence": None,
                "restriction_sites": [],
                "mutations": None,
                "PCR_reactions": {},
                "messages": [],
                "errors": None
            }

            try:
                print(
                    f"Processing sequence {sequence_number}/{len(self.sequencesToDomesticate)}")

                # 1. Preprocess sequence (remove start/stop codons, etc.)
                with debug_context("Preprocessing sequence"):
                    processed_seq, message, _ = self.sequence_preparator.preprocess_sequence(
                        single_seq)

                    if message:
                        sequence_data["messages"].append(message)

                    sequence_data["processed_sequence"] = str(
                        processed_seq) if processed_seq else str(single_seq)

                # 2. Find restriction sites
                with debug_context("Finding restriction sites"):
                    sites_to_mutate = self.sequence_preparator.find_and_summarize_sites(
                        processed_seq, i)
                    sequence_data["restriction_sites"] = sites_to_mutate

                # 3. Mutation analysis and mutation primer design
                mutation_primers = {}

                if sites_to_mutate:
                    with debug_context("Mutation analysis"):
                        mutation_options = self.mutation_analyzer.get_all_mutations(
                            sites_to_mutate)
                        if mutation_options:
                            optimized_mutations, compatibility_matrices = self.mutation_optimizer.optimize_mutations(
                                sequence=processed_seq, mutation_options=mutation_options
                            )
                            sequence_data["mutations"] = {
                                "all_mutation_options": optimized_mutations,
                                "compatibility": compatibility_matrices
                            }
                        else:
                            optimized_mutations, compatibility_matrices = None, None

                    with debug_context("Mutation primer design"):
                        mutation_primers = self.primer_designer.design_mutation_primers(
                            full_sequence=processed_seq,
                            mutation_sets=optimized_mutations,
                            comp_matrices=compatibility_matrices,
                            primer_name=seq_object.get(
                                "primerName", f"Sequence_{i}")
                        )

                # 4. Generate edge primers
                try:
                    mtk_part_left = seq_object.get("mtkPartLeft", "")
                    mtk_part_right = seq_object.get("mtkPartRight", "")
                    print(
                        f"self.mtk_partend_sequences: {self.mtk_partend_sequences}")

                    # Use the extracted MTK parts directly
                    overhang_5_prime = self.utils.get_mtk_partend_sequence(
                        mtk_part_left, "forward", self.kozak)
                    overhang_3_prime = self.utils.get_mtk_partend_sequence(
                        mtk_part_right, "reverse", self.kozak)

                    print(f"overhang_5_prime: {overhang_5_prime}")
                    print(f"overhang_3_prime: {overhang_3_prime}")

                    # Get the primer name from the sequence object if available
                    primer_name = seq_object.get("primerName", f"Sequence_{i}")

                    edge_primers = self.primer_designer.generate_GG_edge_primers(
                        i, processed_seq, overhang_5_prime, overhang_3_prime, primer_name
                    )

                except Exception as e:
                    logger.error(
                        f"Error generating edge primers for sequence {i}: {e}")
                    sequence_data["errors"] = str(e)
                    edge_primers = {}

                # 5. Group primers into PCR reactions
                print(f"Edge primers: {edge_primers}")
                print(f"Mutation primers: {mutation_primers}")
                primers = {**edge_primers, **mutation_primers}
                print(f"Primers: {primers}")
                sequence_data["PCR_reactions"] = self.group_primers_into_pcr_reactions(
                    primers)

                # Store the processed sequence data for this sequence number
                result_data[sequence_number] = sequence_data

            except Exception as e:
                logger.exception(
                    f"Unhandled error processing sequence {i}: {str(e)}")
                sequence_data["errors"] = str(e)
                result_data[sequence_number] = sequence_data

        # Convert any remaining non-serializable objects and return the dictionary
        return self.utils.convert_non_serializable(result_data)

    def _save_primers_to_tsv(self, primer_data: List[List[str]], output_tsv_path: str) -> None:
        """Saves primer data to a TSV file."""
        with debug_context("save_primers_to_tsv"):
            if not primer_data:
                logger.warning("No primer data to save.")
                return

            try:
                with open(output_tsv_path, "w") as tsv_file:
                    tsv_file.write("Primer Name\tSequence\tAmplicon\n")
                    for row in primer_data:
                        tsv_file.write("\t".join(map(str, row)) + "\n")
            except IOError as e:
                logger.error(f"Error writing to file {output_tsv_path}: {e}")
                raise

    def group_primers_into_pcr_reactions(self, primers: Dict[str, List[str]]) -> Dict[str, Dict[str, List[str]]]:
        """
        Organizes primers into PCR reactions following these rules:
        - First mutation reverse primer pairs with edge forward primer (Reaction_1).
        - Each subsequent mutation primer pair forms a new reaction.
        - If no mutation primers exist, only an edge primer reaction is created.

        Args:
            primers (dict): Dictionary containing mutation primers and edge primers.

        Returns:
            dict: PCR reactions grouped properly.
        """
        pcr_reactions = {}
        reaction_counter = 1

        mutation_primers = primers.get("mutation_primers", [])
        edge_primers = primers.get("edge_primers", [])

        if not edge_primers:  # Edge primers are always expected, but handle edge cases
            logger.warning(
                "No edge primers found; PCR reaction assignment may be incomplete.")

        # If no mutations, just assign edge primers to one reaction
        if not mutation_primers:
            pcr_reactions[f"Reaction_{reaction_counter}"] = {
                "edge_primers": edge_primers
            }
            return pcr_reactions

        # First reaction: Mutation reverse primer + Edge forward primer
        pcr_reactions[f"Reaction_{reaction_counter}"] = {
            # First mutation reverse primer
            "mutation_primers": [mutation_primers[0]],
            "edge_primers": [edge_primers[0]]  # Edge forward primer
        }
        reaction_counter += 1

        # Subsequent reactions: Mutation forward + Next mutation reverse
        for i in range(1, len(mutation_primers) - 1, 2):  # Step through pairs
            pcr_reactions[f"Reaction_{reaction_counter}"] = {
                "mutation_primers": [mutation_primers[i], mutation_primers[i + 1]]
            }
            reaction_counter += 1

        return pcr_reactions
