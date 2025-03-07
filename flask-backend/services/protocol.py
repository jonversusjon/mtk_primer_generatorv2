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
from models.primer import Primer, PrimerSet


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
        self.primer_designer = PrimerDesigner(
            kozak=kozak, verbose=verbose, debug=True)
        self.primer_selector = PrimerSelector()
        self.mutation_optimizer = MutationOptimizer(
            verbose=verbose, debug=True)
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

        for idx, seq_object in enumerate(self.sequencesToDomesticate):
            single_seq, mtk_part_left, mtk_part_right, primer_name = getSequenceData(
                seq_object, idx+1)

            sequence_data = {
                "sequence_index": idx,
                "processed_sequence": None,
                "restriction_sites": [],
                "mutations": None,
                "PCR_reactions": {},
                "messages": [],
                "errors": None
            }

            print(
                f"Processing sequence {idx+1}/{len(self.sequencesToDomesticate)}")

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
                    processed_seq, idx)
                sequence_data["restriction_sites"] = sites_to_mutate

            # 3. Mutation analysis and mutation primer design
            mutation_primers = {}

            if sites_to_mutate:

                with debug_context("Mutation analysis"):
                    mutation_options = self.mutation_analyzer.get_all_mutations(
                        sites_to_mutate)
                    optimized_mutations, compatibility_matrices = None, None

                    if mutation_options:
                        optimized_mutations, compatibility_matrices = self.mutation_optimizer.optimize_mutations(
                            sequence=processed_seq, mutation_options=mutation_options
                        )

                        sequence_data["mutations"] = {
                            "all_mutation_options": optimized_mutations,
                            "compatibility": compatibility_matrices
                        }

                with debug_context("Mutation primer design"):
                    mutation_primers = self.primer_designer.design_mutation_primers(
                        full_sequence=processed_seq,
                        mutation_sets=optimized_mutations,
                        comp_matrices=compatibility_matrices,
                        primer_name=seq_object.get(
                            "primerName", f"Sequence_{idx+1}")
                    )
                    sequence_data["mutation_primers"] = mutation_primers

                    # Print designed primers
                    # if mutation_primers:
                    #     print("\n=== DESIGNED PRIMERS ===")
                    #     for i, primer in enumerate(mutation_primers):
                    #         print(
                    #             f"Mutation Primer {i+1} for site {primer.site} at position {primer.position}:")
                    #         print(f"  Forward: {primer.forward.name}")
                    #         print(
                    #             f"  Forward Sequence: {primer.forward.sequence}")
                    #         print(f"  Reverse: {primer.reverse.name}")
                    #         print(
                    #             f"  Reverse Sequence: {primer.reverse.sequence}")
                    #     print("=== END DESIGNED PRIMERS ===\n")

            # 4. Generate edge primers
            edge_primers = self.primer_designer.generate_GG_edge_primers(
                idx, processed_seq, mtk_part_left, mtk_part_right, primer_name
            )
            sequence_data["edge_primers"] = edge_primers

            # Print edge primers
            # if edge_primers:
            #     print("\n=== EDGE PRIMERS ===")
            #     print(f"Forward: {edge_primers['forward_primer']['name']}")
            #     print(
            #         f"Forward Sequence: {edge_primers['forward_primer']['sequence']}")
            #     print(f"Reverse: {edge_primers['reverse_primer']['name']}")
            #     print(
            #         f"Reverse Sequence: {edge_primers['reverse_primer']['sequence']}")
            #     print("=== END EDGE PRIMERS ===\n")

            # 5. Group primers into PCR reactions
            sequence_data["PCR_reactions"] = self.group_primers_into_pcr_reactions(
                sequence_data)

            # Store the processed sequence data for this sequence number
            result_data[idx] = sequence_data

        # Convert any remaining non-serializable objects and return the dictionary
        # for key, value in result_data.items():
        #     for k, v in value.items():
        #         print(f"\n\n{k}: {v}")
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

    def group_primers_into_pcr_reactions(self, sequence_data: Dict) -> Dict[str, Dict[str, Primer]]:
        """
        Groups primers into PCR reactions using chaining logic.

        If no mutation primers exist, a single reaction is created:
        Reaction 1: edge.forward + edge.reverse.

        If mutation primers exist (assumed sorted by position), then:
        Reaction 1: edge.forward + first mutation's reverse
        Reaction 2..n: previous mutation's forward + current mutation's reverse
        Final Reaction: last mutation's forward + edge.reverse
        """
        reactions = {}
        reaction_num = 1

        edge_fw = sequence_data["edge_primers"]["forward_primer"]["sequence"]
        edge_rv = sequence_data["edge_primers"]["reverse_primer"]["sequence"]

        # Extract mutation sets from the sequence data
        mutation_sets = sequence_data.get("mutation_primers", {})

        # Sort mutations by position (ascending order)
        mutations = sorted(mutation_sets, key=lambda m: m.position)

        # Case with no mutations:
        if not mutations:
            reactions[f"Reaction_{reaction_num}"] = {
                "forward": edge_fw, "reverse": edge_rv}
            return reactions

        # Reaction 1: edge forward with the first mutation's reverse primer
        reactions[f"Reaction_{reaction_num}"] = {
            "forward": edge_fw, "reverse": mutations[0].reverse}
        reaction_num += 1

        # Intermediate reactions: chain mutation primers
        for i in range(1, len(mutations)):
            reactions[f"Reaction_{reaction_num}"] = {
                "forward": mutations[i - 1].forward,
                "reverse": mutations[i].reverse
            }
            reaction_num += 1

        # Final Reaction: last mutation's forward primer with edge reverse primer
        reactions[f"Reaction_{reaction_num}"] = {
            "forward": mutations[-1].forward, "reverse": edge_rv}

        return reactions


def getSequenceData(seq_object, i):
    """Extracts sequence data from the provided dictionary."""
    single_seq = seq_object.get("sequence", "")
    mtk_part_left = seq_object.get("mtkPartLeft", "")
    mtk_part_right = seq_object.get("mtkPartRight", "")
    primer_name = seq_object.get("primerName", f"Sequence_{i}")

    return single_seq, mtk_part_left, mtk_part_right, primer_name
