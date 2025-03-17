# services/reactions.py

from typing import Dict
from models import Primer
from log_utils import logger


class ReactionOrganizer():
    """

    """
    
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
        reactions_all = {}
        logger.log_step("Group PCR Reactions", "Starting grouping of primers into PCR reactions.")

        # Retrieve edge primers.
        edge_fw: Primer = sequence_data["edge_primers"]["forward_primer"]
        edge_rv: Primer = sequence_data["edge_primers"]["reverse_primer"]
        logger.log_step("Edge Primers Retrieved", "Edge primers obtained.",
                        {"edge_forward": edge_fw["sequence"], "edge_reverse": edge_rv["sequence"]})

        # Retrieve mutation primers data.
        mutation_primers_data = sequence_data.get("mutation_primers", {})

        # No mutation primers: create default edge-only reaction.
        if not mutation_primers_data:
            default_reaction = {"Reaction_1": {"forward": edge_fw["sequence"], "reverse": edge_rv["sequence"]}}
            reactions_all["default"] = default_reaction
            logger.log_step("No Mutation Primers", "No mutation primers found; created edge-only reaction.", default_reaction)
            return reactions_all

        # Process each mutation set.
        for set_key, solutions in mutation_primers_data.items():
            # Check if the solutions are already nested (a list of solutions) or a flat list.
            if solutions and not isinstance(solutions[0], list):
                logger.log_step("Normalize Solutions", f"Mutation set {set_key} received as flat list; normalizing.")
                solutions = [solutions]
            reactions_all[set_key] = {}
            logger.log_step("Process Mutation Set", f"Processing mutation set {set_key}", {"solution_count": len(solutions)})

            # Process each solution for the current mutation set.
            for sol_index, mutation_solution in enumerate(solutions):
                reaction_dict = {}
                reaction_num = 1

                # If no mutation primers exist for this solution, create a single edge-only reaction.
                if not mutation_solution:
                    # Create a default reaction using edge primers.
                    reaction_label = f"Reaction_{reaction_num}"
                    reaction_dict[reaction_label] = {"forward": edge_fw["sequence"], "reverse": edge_rv["sequence"]}
                    logger.log_step("PCR Reaction Created",
                                    f"{reaction_label} (edge-only) created for set {set_key}, solution {sol_index}",
                                    {"forward": edge_fw["sequence"], "reverse": edge_rv["sequence"]})
                else:
                    # Sort mutation primers by position.
                    mutations = sorted(mutation_solution, key=lambda m: m.position)
                    logger.log_step("Sorted Mutations", f"Sorted mutation primers for set {set_key}, solution {sol_index}",
                                    {"sorted_positions": [m.position for m in mutations]})

                    # Reaction 1: edge forward with the first mutation's reverse primer.
                    reaction_label = f"Reaction_{reaction_num}"
                    reaction_dict[reaction_label] = {"forward": edge_fw["sequence"], "reverse": mutations[0].reverse.sequence}
                    logger.log_step("PCR Reaction Created",
                                    f"{reaction_label} created for set {set_key}, solution {sol_index}",
                                    {"forward": edge_fw["sequence"],
                                     "reverse": mutations[0].reverse.sequence,
                                     "mutation_reverse_position": mutations[0].position})
                    reaction_num += 1

                    # Chain intermediate mutation primers.
                    for i in range(1, len(mutations)):
                        reaction_label = f"Reaction_{reaction_num}"
                        forward_seq = mutations[i - 1].forward.sequence
                        reverse_seq = mutations[i].reverse.sequence
                        reaction_dict[reaction_label] = {"forward": forward_seq, "reverse": reverse_seq}
                        logger.log_step("PCR Reaction Created",
                                        f"{reaction_label} created for set {set_key}, solution {sol_index}",
                                        {"forward": forward_seq,
                                         "reverse": reverse_seq,
                                         "from_mutation_position": mutations[i - 1].position,
                                         "to_mutation_position": mutations[i].position})
                        reaction_num += 1

                    # Final Reaction: last mutation's forward with edge reverse.
                    reaction_label = f"Reaction_{reaction_num}"
                    final_forward = mutations[-1].forward.sequence
                    reaction_dict[reaction_label] = {"forward": final_forward, "reverse": edge_rv["sequence"]}
                    logger.log_step("PCR Reaction Created",
                                    f"{reaction_label} created for set {set_key}, solution {sol_index}",
                                    {"forward": final_forward,
                                     "reverse": edge_rv["sequence"],
                                     "from_mutation_position": mutations[-1].position})

                # Save the PCR reaction grouping for this solution.
                reactions_all[set_key][sol_index] = reaction_dict
                logger.log_step("Mutation Set Processed",
                                f"Completed grouping for set {set_key}, solution {sol_index}",
                                {"total_reactions": reaction_num})

        logger.log_step("Group PCR Reactions Complete", "Completed grouping of all PCR reactions.",
                        {"total_mutation_sets": len(reactions_all)})
        return reactions_all
