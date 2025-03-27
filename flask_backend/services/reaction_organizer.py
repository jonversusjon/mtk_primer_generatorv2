# services/reactions.py

from typing import List

from flask_backend.logging import logger
from flask_backend.models import DomesticationResult, Primer, MutationPrimerSet, PCRReaction

class ReactionOrganizer():
    """
    Organizes mutation primer sets into PCR reactions.
    """

    def calculate_amplicon_size(self, forward_primer: Primer, reverse_primer: Primer) -> int:
        """
        Calculate the amplicon size based on primer positions.
        This is a placeholder - replace with actual implementation.
        """
        # In a real implementation, you would calculate this based on 
        # the positions of the primers on the template
        return 500  # Placeholder value

    def create_pcr_reaction(self, name: str, forward_primer: Primer, reverse_primer: Primer) -> PCRReaction:
        """
        Create a PCRReaction object with all required fields.
        """
        amplicon_size = self.calculate_amplicon_size(forward_primer, reverse_primer)
        
        return PCRReaction(
            name=name,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            amplicon_size=amplicon_size
        )

    def group_primers_into_pcr_reactions(
        self,
        domestication_result: DomesticationResult,
        send_update: callable) -> dict:
        """
        Groups primers into nested PCR reactions using chaining logic for each mutation solution.

        Expects domestication_result to have:
        - edge_primers: an object with properties "forward" and "reverse" (both Primer objects)
        - mut_primers: a list of MutationPrimerSet objects,
                        each containing mut_primer_pairs (list of MutationPrimerPair objects).

        The nested dictionary structure returned is:
            {
                "mutation_sets": [
                    {
                        "set_id": int,
                        "solutions": [
                            {
                                "solution_id": int,
                                "reactions": [
                                    {
                                        "reaction_id": str,    # same as reaction name
                                        "forward_primer": <Primer>,
                                        "reverse_primer": <Primer>
                                    },
                                    ...
                                ]
                            }
                        ]
                    },
                    ...
                ]
            }

        For each mutation set and solution:
        - If no mutation primers exist, a single edge-only reaction is created.
        - If mutation primers exist (assumed sorted by position), then:
                Reaction 1: edge forward + first mutation's reverse primer
                Reaction 2..n: previous mutation's forward primer + current mutation's reverse primer
                Final Reaction: last mutation's forward primer + edge reverse

        Returns:
            A nested dictionary representing all reactions across all mutation sets.
        """
        nested_reactions = {"mutation_sets": []}
        logger.log_step("Group PCR Reactions", "Starting grouping of primers into nested PCR reactions.")

        # Retrieve edge primers.
        edge_fw: Primer = domestication_result.edge_primers.forward
        edge_rv: Primer = domestication_result.edge_primers.reverse
        logger.log_step("Edge Primers Retrieved", "Edge primers obtained.",
                        {"edge_forward": edge_fw.sequence, "edge_reverse": edge_rv.sequence})

        # Retrieve mutation primers data.
        mut_primers_collection: List[MutationPrimerSet] = domestication_result.mut_primers
        logger.log_step("Mutation Primers Retrieved", "Mutation primers data obtained.",
                        {"mutation_set_count": len(mut_primers_collection)})
        logger.log_step("Mutation Primers Data", "Mutation primers data details.",
                        {"mutation_primers": mut_primers_collection})

        # If no mutation primers, create a default mutation set/solution with an edge-only reaction.
        if not mut_primers_collection:
            reaction_label = "reaction_1"
            pcr_reaction = self.create_pcr_reaction(
                name=reaction_label,
                forward_primer=edge_fw,
                reverse_primer=edge_rv
            )
            logger.log_step("No Mutation Primers",
                            "No mutation primers found; created edge-only reaction.",
                            {"reaction": reaction_label})
            default_mutation_set = {
                "set_id": 0,
                "solutions": [
                    {
                        "solution_id": 0,
                        "reactions": [pcr_reaction]
                    }
                ]
            }
            nested_reactions["mutation_sets"].append(default_mutation_set)
            return nested_reactions

        # Process each mutation set.
        for set_idx, mut_primer_set in enumerate(mut_primers_collection):
            logger.log_step("Process Mutation Set", f"Processing mutation set {set_idx}/{len(mut_primers_collection)}", "")
            logger.log_step("Mutation Primer Set", "Mutation primer set details.",
                            {"mutation_primer_set": mut_primer_set})

            # Initialize a solution object (currently one per mutation set).
            solution_data = {"solution_id": 0, "reactions": []}
            reaction_num = 1

            mut_primer_pairs = mut_primer_set.mut_primer_pairs

            # If no mutation primers exist for this solution, create a single edge-only reaction.
            if not mut_primer_pairs:
                reaction_label = f"set{set_idx}_sol0_reaction_{reaction_num}"
                pcr_reaction = self.create_pcr_reaction(
                    name=reaction_label,
                    forward_primer=edge_fw,
                    reverse_primer=edge_rv
                )
                solution_data["reactions"].append(pcr_reaction)
                logger.log_step("PCR Reaction Created",
                                f"{reaction_label} (edge-only) created for set {set_idx}, solution 0",
                                {"forward": edge_fw.sequence, "reverse": edge_rv.sequence})
            else:
                logger.log_step("Mutation Primer Pairs", "Data available",
                                f"mut_primer_pairs count: {len(mut_primer_pairs)} pairs")
                # Sort mutation primers by position.
                ordered_mut_primer_pairs = sorted(mut_primer_pairs, key=lambda k: k.position)

                # Reaction 1: edge forward with the first mutation's reverse primer.
                reaction_label = f"set{set_idx}_sol0_reaction_{reaction_num}"
                first_rev_primer = ordered_mut_primer_pairs[0].reverse
                pcr_reaction = self.create_pcr_reaction(
                    name=reaction_label,
                    forward_primer=edge_fw,
                    reverse_primer=first_rev_primer
                )
                solution_data["reactions"].append(pcr_reaction)
                logger.log_step("PCR Reaction Created",
                                f"{reaction_label} created for set {set_idx}, solution 0",
                                {"forward": edge_fw.sequence,
                                "reverse": first_rev_primer.sequence,
                                "mutation_reverse_position": ordered_mut_primer_pairs[0].position})
                reaction_num += 1

                # Chain intermediate mutation primers.
                for i in range(1, len(ordered_mut_primer_pairs)):
                    reaction_label = f"set{set_idx}_sol0_reaction_{reaction_num}"
                    forward_primer = ordered_mut_primer_pairs[i - 1].forward
                    reverse_primer = ordered_mut_primer_pairs[i].reverse
                    pcr_reaction = self.create_pcr_reaction(
                        name=reaction_label,
                        forward_primer=forward_primer,
                        reverse_primer=reverse_primer
                    )
                    solution_data["reactions"].append(pcr_reaction)
                    logger.log_step("PCR Reaction Created",
                                    f"{reaction_label} created for set {set_idx}, solution 0",
                                    {"forward": forward_primer.sequence,
                                    "reverse": reverse_primer.sequence,
                                    "from_mutation_position": ordered_mut_primer_pairs[i - 1].position,
                                    "to_mutation_position": ordered_mut_primer_pairs[i].position})
                    reaction_num += 1

                # Final Reaction: last mutation's forward with edge reverse.
                reaction_label = f"set{set_idx}_sol0_reaction_{reaction_num}"
                last_forward_primer = ordered_mut_primer_pairs[-1].forward
                pcr_reaction = self.create_pcr_reaction(
                    name=reaction_label,
                    forward_primer=last_forward_primer,
                    reverse_primer=edge_rv
                )
                solution_data["reactions"].append(pcr_reaction)
                logger.log_step("PCR Reaction Created",
                                f"{reaction_label} created for set {set_idx}, solution 0",
                                {"forward": last_forward_primer.sequence,
                                "reverse": edge_rv.sequence,
                                "from_mutation_position": ordered_mut_primer_pairs[-1].position})

            logger.log_step("Mutation Set Processed",
                            f"Completed grouping for set {set_idx}, solution 0",
                            {"total_reactions": reaction_num})

            mutation_set_data = {
                "set_id": set_idx,
                "solutions": [solution_data]
            }
            nested_reactions["mutation_sets"].append(mutation_set_data)

        logger.log_step("Group PCR Reactions Complete", "Completed grouping of all nested PCR reactions.",
                        {"total_mutation_sets": len(nested_reactions["mutation_sets"])})
        send_update(message="PCR Reaction Grouping Complete", prog=100)
        
        return nested_reactions