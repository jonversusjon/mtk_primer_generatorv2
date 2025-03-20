# services/reactions.py

from typing import Dict, List
from models import Primer
from log_utils import logger
from models import DomesticationResult, MutationPrimerSet

class ReactionOrganizer():
    """
    Organizes mutation primer sets into PCR reactions.
    """

    def group_primers_into_pcr_reactions(self, domestication_result: DomesticationResult) -> Dict:
        """
        Groups primers into PCR reactions using chaining logic for each mutation solution.
        
        Expects domestication_result to have:
          - edge_primers: an object with properties "forward" and "reverse" (both Primer objects)
          - mut_primers: a dict mapping mutation set indices to lists of solutions,
                         where each solution is a list of MutationPrimer objects.
        
        For each mutation set and solution:
         - If no mutation primers exist, creates a single reaction:
                Reaction_1: edge.forward + edge.reverse.
         - If mutation primers exist (assumed sorted by position), then:
                Reaction_1: edge forward + first mutation's reverse primer
                Reaction_2..n: previous mutation's forward primer + current mutation's reverse primer
                Final Reaction: last mutation's forward primer + edge reverse
        
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
        edge_fw: Primer = domestication_result.edge_primers.forward
        edge_rv: Primer = domestication_result.edge_primers.reverse
        logger.log_step("Edge Primers Retrieved", "Edge primers obtained.",
                        {"edge_forward": edge_fw.sequence, "edge_reverse": edge_rv.sequence})

        # Retrieve mutation primers data.
        mut_primers_collection: List[MutationPrimerSet] = domestication_result.mut_primers
        logger.log_step("Mutation Primers Retrieved", "Mutation primers data obtained.",
                        {" ****** mutation_set ****** ": len(mut_primers_collection)})
        logger.log_step("Mutation Primers Data", "Mutation primers data details.",
                        {"mutation_primers": mut_primers_collection})
        
        # No mutation primers: create default edge-only reaction.
        if not mut_primers_collection:
            default_reaction = {"Reaction_1": {"forward": edge_fw.sequence, "reverse": edge_rv.sequence}}
            reactions_all["default"] = default_reaction
            logger.log_step("No Mutation Primers", "No mutation primers found; created edge-only reaction.", default_reaction)
            return reactions_all
        
        # Process each mutation set.
        for idx, mut_primer_set in enumerate(mut_primers_collection):
            reactions_all[idx] = {}
            logger.log_step("Process Mutation Set", f"Processing mutation set {idx}/{len(mut_primers_collection)}", "")
            logger.log_step("Mutation Primer Set", "Mutation primer set details.",
                            {"mutation_primer_set": mut_primer_set})
            
            # Get mutation primer pairs from the MutationPrimerSet
            mut_primer_pairs = mut_primer_set.mut_primer_pairs
            
            # Process this solution
            mut_primer_idx = 0  # Single solution per mutation set
            reaction_dict = {}
            reaction_num = 1

            # If no mutation primers exist for this solution, create a single edge-only reaction.
            if not mut_primer_pairs:
                reaction_label = f"Reaction_{reaction_num}"
                reaction_dict[reaction_label] = {"forward": edge_fw.sequence, "reverse": edge_rv.sequence}
                logger.log_step("PCR Reaction Created",
                            f"{reaction_label} (edge-only) created for set {idx}, solution {mut_primer_idx}",
                            {"forward": edge_fw.sequence, "reverse": edge_rv.sequence})
            else:
                logger.log_step("Mutation primer pairs", "data", f"mut_primer_pairs: {len(mut_primer_pairs)} pairs")
                
                # Sort mutation primers by position
                ordered_mut_primer_pairs = sorted(mut_primer_pairs, key=lambda k: k.position)
                
                # Reaction 1: edge forward with the first mutation's reverse primer.
                reaction_label = f"Reaction_{reaction_num}"
                reaction_dict[reaction_label] = {"forward": edge_fw.sequence, "reverse": ordered_mut_primer_pairs[0].reverse.sequence}
                logger.log_step("PCR Reaction Created",
                            f"{reaction_label} created for set {idx}, solution {mut_primer_idx}",
                            {"forward": edge_fw.sequence,
                                "reverse": ordered_mut_primer_pairs[0].reverse.sequence,
                                "mutation_reverse_position": ordered_mut_primer_pairs[0].position})
                reaction_num += 1

                # Chain intermediate mutation primers.
                for i in range(1, len(ordered_mut_primer_pairs)):
                    reaction_label = f"Reaction_{reaction_num}"
                    forward_seq = ordered_mut_primer_pairs[i - 1].forward.sequence
                    reverse_seq = ordered_mut_primer_pairs[i].reverse.sequence
                    reaction_dict[reaction_label] = {"forward": forward_seq, "reverse": reverse_seq}
                    logger.log_step("PCR Reaction Created",
                                    f"{reaction_label} created for set {idx}, solution {mut_primer_idx}",
                                    {"forward": forward_seq,
                                    "reverse": reverse_seq,
                                    "from_mutation_position": ordered_mut_primer_pairs[i - 1].position,
                                    "to_mutation_position": ordered_mut_primer_pairs[i].position})
                    reaction_num += 1

                # Final Reaction: last mutation's forward with edge reverse.
                reaction_label = f"Reaction_{reaction_num}"
                final_forward = ordered_mut_primer_pairs[-1].forward.sequence
                reaction_dict[reaction_label] = {"forward": final_forward, "reverse": edge_rv.sequence}
                logger.log_step("PCR Reaction Created",
                                f"{reaction_label} created for set {idx}, solution {mut_primer_idx}",
                                {"forward": final_forward,
                                "reverse": edge_rv.sequence,
                                "from_mutation_position": ordered_mut_primer_pairs[-1].position})

            # Save the PCR reaction grouping for this solution.
            reactions_all[idx][mut_primer_idx] = reaction_dict
            logger.log_step("Mutation Set Processed",
                            f"Completed grouping for set {idx}, solution {mut_primer_idx}",
                            {"total_reactions": reaction_num})

        logger.log_step("Group PCR Reactions Complete", "Completed grouping of all PCR reactions.",
                        {"total_mutation_sets": len(reactions_all)})
        
        return reactions_all

