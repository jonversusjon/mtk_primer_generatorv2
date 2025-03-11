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
from models.primer import Primer, PrimerSet, MutationPrimer
from services.debug.debug_utils import MutationDebugger, visualize_matrix
from services.debug.debug_mixin import DebugMixin


class GoldenGateProtocol(DebugMixin):
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
        max_results: int = 1,
        verbose: bool = False
    ):
        self.logger = logger.getChild("GoldenGateProtocol")

        self.utils = GoldenGateUtils()
        self.sequence_preparator = SequencePreparator()
        self.mutation_analyzer = MutationAnalyzer(
            codon_usage_dict=codon_usage_dict,
            max_mutations=max_mutations,
            verbose=verbose,
            debug=False,
        )
        self.mutation_optimizer = MutationOptimizer(
            verbose=verbose, debug=False)
        self.primer_designer = PrimerDesigner(
            kozak=kozak, verbose=verbose, debug=True)

        self.logger.debug(
            f"GoldenGateProtocol initialized with codon_usage_dict: {codon_usage_dict}")
        if verbose:
            self.logger.info("GoldenGateProtocol is running in verbose mode.")

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

        self.max_results = max_results

    def create_gg_protocol(self) -> dict:
        """
        Main function to orchestrate the Golden Gate protocol creation.
        Returns:
            dict: A dictionary containing protocol details.
        """
        print(f"Starting Golden Gate protocol creation...")
        logger.info("Starting Golden Gate protocol creation...")

        result_data = {}

        for idx, seq_object in enumerate(self.sequencesToDomesticate):
            single_seq, mtk_part_left, mtk_part_right, primer_name = getSequenceData(
                seq_object, idx+1)

            sequence_data = {
                "sequence_index": idx,
                "processed_sequence": None,
                "mtk_part_left": mtk_part_left,
                "mtk_part_right": mtk_part_right,
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
                    single_seq, mtk_part_left)

                if message:
                    sequence_data["messages"].append(message)

                sequence_data["processed_sequence"] = str(
                    processed_seq) if processed_seq else str(single_seq)

            # 2. Find restriction sites
            with debug_context("Finding restriction sites"):
                sites_to_mutate = self.sequence_preparator.find_sites_to_mutate(
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
                            mutation_options=mutation_options
                        )

                        # TODO: add valid solution tracking to enable retuning multiple solutions up to
                        # max_results number of solutions to the user
                        sequence_data["mutations"] = {
                            "all_mutation_options": optimized_mutations,
                            "compatibility": compatibility_matrices
                        }

                with debug_context("Mutation primer design"):

                    mutation_primers = self.primer_designer.design_mutation_primers(
                        mutation_sets=optimized_mutations,
                        comp_matrices=compatibility_matrices,
                        primer_name=seq_object.get(
                            "primerName", f"Sequence_{idx+1}")
                    )
                    sequence_data["mutation_primers"] = mutation_primers

            # 4. Generate edge primers
            edge_primers = self.primer_designer.generate_GG_edge_primers(
                idx, processed_seq, mtk_part_left, mtk_part_right, primer_name
            )
            sequence_data["edge_primers"] = edge_primers

            # 5. Group primers into PCR reactions
            sequence_data["PCR_reactions"] = self.group_primers_into_pcr_reactions(
                sequence_data)

            # Store the processed sequence data for this sequence number
            result_data[idx] = sequence_data

        return self.utils.convert_non_serializable(result_data)

    def group_primers_into_pcr_reactions(self, sequence_data: Dict) -> Dict[str, Dict[str, Primer]]:
        """
        Groups primers into PCR reactions using chaining logic.

        If no mutation primers exist, a single reaction is created:
        Reaction 1: edge.forward + edge.reverse.

        If mutation primers exist (assumed sorted by position), then:
        Reaction 1: edge.forward + first mutation's reverse primer
        Reaction 2..n: previous mutation's forward primer + current mutation's reverse primer
        Final Reaction: last mutation's forward primer + edge.reverse
        """
        reactions = {}
        reaction_num = 1

        self.log_step("Group PCR Reactions",
                      "Starting grouping of primers into PCR reactions.")

        # Retrieve edge primers: Adjusted to match expected structure.
        edge_fw: Primer = sequence_data["edge_primers"]["forward_primer"]
        edge_rv: Primer = sequence_data["edge_primers"]["reverse_primer"]

        self.log_step("Edge Primers Retrieved", "Edge primers obtained.",
                      {"edge_forward": edge_fw["sequence"], "edge_reverse": edge_rv["sequence"]})

        # Get the mutation primers list (assumed to be a list of MutationPrimer instances).
        mutation_primers: list[MutationPrimer] = sequence_data.get(
            "mutation_primers", [])
        self.log_step("Mutation Primers Retrieved", "Retrieved mutation primers.",
                      {"mutation_primer_count": len(mutation_primers)})

        # If no mutation primers exist, create a single reaction with just the edge primers.
        if not mutation_primers:
            self.log_step("No Mutation Primers",
                          "No mutation primers found; creating single edge-only reaction.")
            reactions[f"Reaction_{reaction_num}"] = {
                "forward": edge_fw["sequence"], "reverse": edge_rv["sequence"]}
            self.log_step("PCR Reaction Created", f"Reaction_{reaction_num} created.",
                          {"forward": edge_fw["sequence"], "reverse": edge_rv["sequence"]})
            return reactions

        # Sort mutation primers by position.
        mutations = sorted(mutation_primers, key=lambda m: m.position)
        self.log_step("Sorted Mutation Primers", "Sorted mutation primers by position.",
                      {"sorted_positions": [m.position for m in mutations]})

        # Reaction 1: edge forward with first mutation's reverse primer.
        reaction_label = f"Reaction_{reaction_num}"
        reactions[reaction_label] = {
            "forward": edge_fw["sequence"], "reverse": mutations[0].reverse.sequence}
        self.log_step("PCR Reaction Created", f"{reaction_label} created.",
                      {"forward": edge_fw["sequence"], "reverse": mutations[0].reverse.sequence,
                       "edge_forward": True, "mutation_reverse_position": mutations[0].position})
        reaction_num += 1

        # Chain intermediate mutation primers.
        for i in range(1, len(mutations)):
            reaction_label = f"Reaction_{reaction_num}"
            forward_seq = mutations[i - 1].forward.sequence
            reverse_seq = mutations[i].reverse.sequence
            reactions[reaction_label] = {
                "forward": forward_seq, "reverse": reverse_seq}
            self.log_step("PCR Reaction Created", f"{reaction_label} created.",
                          {"forward": forward_seq, "reverse": reverse_seq,
                           "from_mutation_position": mutations[i - 1].position,
                           "to_mutation_position": mutations[i].position})
            reaction_num += 1

        # Final Reaction: last mutation's forward with edge reverse.
        reaction_label = f"Reaction_{reaction_num}"
        final_forward = mutations[-1].forward.sequence
        reactions[reaction_label] = {
            "forward": final_forward, "reverse": edge_rv["sequence"]}
        self.log_step("PCR Reaction Created", f"{reaction_label} created.",
                      {"forward": final_forward, "reverse": edge_rv["sequence"],
                       "from_mutation_position": mutations[-1].position, "edge_reverse": True})

        self.log_step("Group PCR Reactions Complete", "Completed grouping of PCR reactions.",
                      {"total_reactions": reaction_num})
        return reactions


def getSequenceData(seq_object, i):
    """Extracts sequence data from the provided dictionary."""
    single_seq = seq_object.get("sequence", "")
    mtk_part_left = seq_object.get("mtkPartLeft", "")
    mtk_part_right = seq_object.get("mtkPartRight", "")
    primer_name = seq_object.get("primerName", f"Sequence_{i}")

    return single_seq, mtk_part_left, mtk_part_right, primer_name
