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

    @DebugMixin.debug_wrapper
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

                        sequence_data["mutations"] = {
                            "all_mutation_options": optimized_mutations,
                            "compatibility": compatibility_matrices
                        }

                with debug_context("Mutation primer design"):
                    # TODO: add valid solution tracking to enable returning multiple solutions up to
                    # max_results number of solutions to the user
                    mutation_primers = self.primer_designer.design_mutation_primers(
                        mutation_sets=optimized_mutations,
                        comp_matrices=compatibility_matrices,
                        primer_name=seq_object.get(
                            "primerName", f"Sequence_{idx+1}",
                        ),
                        max_results=self.max_results,
                    )
                    sequence_data["mutation_primers"] = mutation_primers

            print(f"Mutation primers designed: {mutation_primers}")
            print(f"designing edge primers...")
            # 4. Generate edge primers
            edge_primers = self.primer_designer.generate_GG_edge_primers(
                idx, processed_seq, mtk_part_left, mtk_part_right, primer_name
            )
            sequence_data["edge_primers"] = edge_primers

            # 5. Group primers into PCR reactions
            print(f"Grouping primers into PCR reactions...")
            self.log_step("Group PCR Reactions",
                          "sequence_data: {sequence_data}")
            sequence_data["PCR_reactions"] = self.group_primers_into_pcr_reactions(
                sequence_data)
            print(f"Finished grouping primers into PCR reactions...")
            # Store the processed sequence data for this sequence number
            result_data[idx] = sequence_data

        return self.utils.convert_non_serializable(result_data)

    @DebugMixin.debug_wrapper
    def group_primers_into_pcr_reactions(self, sequence_data: Dict) -> Dict:
        """
        Groups primers into PCR reactions using chaining logic for each mutation solution.

        Expects sequence_data to have:
            "edge_primers": dict with keys "forward_primer" and "reverse_primer"
            "mutation_primers": dict mapping mutation set indices to either:
                                - a list of MutationPrimer objects (for a single solution) or
                                - a list of solutions, where each solution is a list of MutationPrimer objects.

        For each mutation set and solution:
        - If no mutation primers exist, creates a single reaction:
                Reaction_1: edge.forward + edge.reverse.
        - If mutation primers exist (assumed sorted by position), then:
                Reaction_1: edge.forward + first mutation's reverse primer
                Reaction_2..n: previous mutation's forward primer + current mutation's reverse primer
                Final Reaction: last mutation's forward primer + edge.reverse

        Returns:
            A nested dictionary in the form:
            {
                mutation_set_index: {
                    solution_index: {
                        "Reaction_1": {"forward": <seq>, "reverse": <seq>},
                        "Reaction_2": { ... },
                        ...
                    },
                    ...
                },
                ...
            }
        """
        reactions_all = {}  # To hold all PCR reaction groups by mutation set.
        self.log_step("Group PCR Reactions",
                      "Starting grouping of primers into PCR reactions.")

        # Retrieve edge primers.
        edge_fw: Primer = sequence_data["edge_primers"]["forward_primer"]
        edge_rv: Primer = sequence_data["edge_primers"]["reverse_primer"]
        self.log_step("Edge Primers Retrieved", "Edge primers obtained.",
                      {"edge_forward": edge_fw["sequence"], "edge_reverse": edge_rv["sequence"]})

        # Retrieve the mutation primers data.
        mutation_primers_data = sequence_data.get("mutation_primers", {})

        # If there are no mutation primers at all, create a default reaction.
        if not mutation_primers_data:
            reactions_all["default"] = {
                "Reaction_1": {"forward": edge_fw["sequence"], "reverse": edge_rv["sequence"]}
            }
            self.log_step("No Mutation Primers",
                          "No mutation primers found; creating single edge-only reaction.",
                          {"Reaction_1": {"forward": edge_fw["sequence"], "reverse": edge_rv["sequence"]}})
            return reactions_all

        # Process each mutation set.
        for set_key, solutions in mutation_primers_data.items():
            # Check if the solutions are already nested (a list of solutions) or a flat list.
            if solutions and not isinstance(solutions[0], list):
                self.log_step("Group PCR Reactions",
                              f"Mutation set {set_key} received as flat list; wrapping in a list for uniform processing.")
                solutions = [solutions]
            reactions_all[set_key] = {}
            self.log_step("Group PCR Reactions",
                          f"Processing mutation set {set_key}",
                          {"solution_count": len(solutions)})

            # Process each solution for the current mutation set.
            for sol_index, mutation_solution in enumerate(solutions):
                reaction_dict = {}
                reaction_num = 1

                # If no mutation primers exist for this solution, create a single edge-only reaction.
                if not mutation_solution:
                    reaction_label = f"Reaction_{reaction_num}"
                    reaction_dict[reaction_label] = {
                        "forward": edge_fw["sequence"],
                        "reverse": edge_rv["sequence"]
                    }
                    self.log_step("PCR Reaction Created",
                                  f"{reaction_label} created for mutation set {set_key}, solution {sol_index}",
                                  {"forward": edge_fw["sequence"], "reverse": edge_rv["sequence"]})
                else:
                    # Sort the mutation primers by position.
                    mutations = sorted(mutation_solution,
                                       key=lambda m: m.position)
                    self.log_step("Sorted Mutation Primers",
                                  f"Sorted mutation primers for mutation set {set_key}, solution {sol_index}",
                                  {"sorted_positions": [m.position for m in mutations]})

                    # Reaction 1: edge forward with the first mutation's reverse primer.
                    reaction_label = f"Reaction_{reaction_num}"
                    reaction_dict[reaction_label] = {
                        "forward": edge_fw["sequence"],
                        "reverse": mutations[0].reverse.sequence
                    }
                    self.log_step("PCR Reaction Created",
                                  f"{reaction_label} created for mutation set {set_key}, solution {sol_index}",
                                  {"forward": edge_fw["sequence"],
                                   "reverse": mutations[0].reverse.sequence,
                                   "edge_forward": True,
                                   "mutation_reverse_position": mutations[0].position})
                    reaction_num += 1

                    # Chain intermediate mutation primers.
                    for i in range(1, len(mutations)):
                        reaction_label = f"Reaction_{reaction_num}"
                        forward_seq = mutations[i - 1].forward.sequence
                        reverse_seq = mutations[i].reverse.sequence
                        reaction_dict[reaction_label] = {
                            "forward": forward_seq,
                            "reverse": reverse_seq
                        }
                        self.log_step("PCR Reaction Created",
                                      f"{reaction_label} created for mutation set {set_key}, solution {sol_index}",
                                      {"forward": forward_seq,
                                       "reverse": reverse_seq,
                                       "from_mutation_position": mutations[i - 1].position,
                                       "to_mutation_position": mutations[i].position})
                        reaction_num += 1

                    # Final Reaction: last mutation's forward with edge reverse.
                    reaction_label = f"Reaction_{reaction_num}"
                    final_forward = mutations[-1].forward.sequence
                    reaction_dict[reaction_label] = {
                        "forward": final_forward,
                        "reverse": edge_rv["sequence"]
                    }
                    self.log_step("PCR Reaction Created",
                                  f"{reaction_label} created for mutation set {set_key}, solution {sol_index}",
                                  {"forward": final_forward,
                                   "reverse": edge_rv["sequence"],
                                   "from_mutation_position": mutations[-1].position,
                                   "edge_reverse": True})

                # Save the PCR reaction grouping for this solution.
                reactions_all[set_key][sol_index] = reaction_dict
                self.log_step("Group PCR Reactions",
                              f"Completed grouping for mutation set {set_key}, solution {sol_index}",
                              {"total_reactions": reaction_num})

        self.log_step("Group PCR Reactions Complete",
                      "Completed grouping of all PCR reactions.",
                      {"total_mutation_sets": len(reactions_all)})
        return reactions_all


def getSequenceData(seq_object, i):
    """Extracts sequence data from the provided dictionary."""
    single_seq = seq_object.get("sequence", "")
    mtk_part_left = seq_object.get("mtkPartLeft", "")
    mtk_part_right = seq_object.get("mtkPartRight", "")
    primer_name = seq_object.get("primerName", f"Sequence_{i}")

    return single_seq, mtk_part_left, mtk_part_right, primer_name
