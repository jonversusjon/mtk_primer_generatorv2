# services/reactions.py

from typing import List
from log_utils import logger
from models import DomesticationResult, Primer, MutationPrimerSet, PCRReaction

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

    def group_primers_into_pcr_reactions(self, domestication_result: DomesticationResult) -> List[PCRReaction]:
        """
        Groups primers into PCR reactions using chaining logic for each mutation solution.
        
        Expects domestication_result to have:
          - edge_primers: an object with properties "forward" and "reverse" (both Primer objects)
          - mut_primers: a list of MutationPrimerSet objects,
                         each containing mut_primer_pairs (list of MutationPrimerPair objects).
        
        For each mutation set and solution:
         - If no mutation primers exist, creates a single reaction:
                Reaction_1: edge.forward + edge.reverse.
         - If mutation primers exist (assumed sorted by position), then:
                Reaction_1: edge forward + first mutation's reverse primer
                Reaction_2..n: previous mutation's forward primer + current mutation's reverse primer
                Final Reaction: last mutation's forward primer + edge reverse
        
        Returns:
            A list of PCRReaction objects representing all reactions across all mutation sets.
        """
        all_reactions = []
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
            reaction_name = "Reaction_1"
            pcr_reaction = self.create_pcr_reaction(
                name=reaction_name,
                forward_primer=edge_fw,
                reverse_primer=edge_rv
            )
            all_reactions.append(pcr_reaction)
            logger.log_step("No Mutation Primers", "No mutation primers found; created edge-only reaction.", pcr_reaction)
            return all_reactions
        
        # Process each mutation set.
        for idx, mut_primer_set in enumerate(mut_primers_collection):
            logger.log_step("Process Mutation Set", f"Processing mutation set {idx}/{len(mut_primers_collection)}", "")
            logger.log_step("Mutation Primer Set", "Mutation primer set details.",
                            {"mutation_primer_set": mut_primer_set})
            
            # Get mutation primer pairs from the MutationPrimerSet
            mut_primer_pairs = mut_primer_set.mut_primer_pairs
            
            # Process this solution
            mut_primer_idx = 0  # Single solution per mutation set
            reaction_num = 1

            # If no mutation primers exist for this solution, create a single edge-only reaction.
            if not mut_primer_pairs:
                reaction_label = f"Set{idx}_Sol{mut_primer_idx}_Reaction_{reaction_num}"
                pcr_reaction = self.create_pcr_reaction(
                    name=reaction_label,
                    forward_primer=edge_fw,
                    reverse_primer=edge_rv
                )
                all_reactions.append(pcr_reaction)
                logger.log_step("PCR Reaction Created",
                            f"{reaction_label} (edge-only) created for set {idx}, solution {mut_primer_idx}",
                            {"forward": edge_fw.sequence, "reverse": edge_rv.sequence})
            else:
                logger.log_step("Mutation primer pairs", "data", f"mut_primer_pairs: {len(mut_primer_pairs)} pairs")
                
                # Sort mutation primers by position
                ordered_mut_primer_pairs = sorted(mut_primer_pairs, key=lambda k: k.position)
                
                # Reaction 1: edge forward with the first mutation's reverse primer.
                reaction_label = f"Set{idx}_Sol{mut_primer_idx}_Reaction_{reaction_num}"
                first_rev_primer = ordered_mut_primer_pairs[0].reverse
                pcr_reaction = self.create_pcr_reaction(
                    name=reaction_label,
                    forward_primer=edge_fw,
                    reverse_primer=first_rev_primer
                )
                all_reactions.append(pcr_reaction)
                logger.log_step("PCR Reaction Created",
                            f"{reaction_label} created for set {idx}, solution {mut_primer_idx}",
                            {"forward": edge_fw.sequence,
                                "reverse": first_rev_primer.sequence,
                                "mutation_reverse_position": ordered_mut_primer_pairs[0].position})
                reaction_num += 1

                # Chain intermediate mutation primers.
                for i in range(1, len(ordered_mut_primer_pairs)):
                    reaction_label = f"Set{idx}_Sol{mut_primer_idx}_Reaction_{reaction_num}"
                    forward_primer = ordered_mut_primer_pairs[i - 1].forward
                    reverse_primer = ordered_mut_primer_pairs[i].reverse
                    pcr_reaction = self.create_pcr_reaction(
                        name=reaction_label,
                        forward_primer=forward_primer,
                        reverse_primer=reverse_primer
                    )
                    all_reactions.append(pcr_reaction)
                    logger.log_step("PCR Reaction Created",
                                    f"{reaction_label} created for set {idx}, solution {mut_primer_idx}",
                                    {"forward": forward_primer.sequence,
                                    "reverse": reverse_primer.sequence,
                                    "from_mutation_position": ordered_mut_primer_pairs[i - 1].position,
                                    "to_mutation_position": ordered_mut_primer_pairs[i].position})
                    reaction_num += 1

                # Final Reaction: last mutation's forward with edge reverse.
                reaction_label = f"Set{idx}_Sol{mut_primer_idx}_Reaction_{reaction_num}"
                last_forward_primer = ordered_mut_primer_pairs[-1].forward
                pcr_reaction = self.create_pcr_reaction(
                    name=reaction_label,
                    forward_primer=last_forward_primer,
                    reverse_primer=edge_rv
                )
                all_reactions.append(pcr_reaction)
                logger.log_step("PCR Reaction Created",
                                f"{reaction_label} created for set {idx}, solution {mut_primer_idx}",
                                {"forward": last_forward_primer.sequence,
                                "reverse": edge_rv.sequence,
                                "from_mutation_position": ordered_mut_primer_pairs[-1].position})

            logger.log_step("Mutation Set Processed",
                            f"Completed grouping for set {idx}, solution {mut_primer_idx}",
                            {"total_reactions": reaction_num})

        logger.log_step("Group PCR Reactions Complete", "Completed grouping of all PCR reactions.",
                        {"total_reactions": len(all_reactions)})
        
        return all_reactions