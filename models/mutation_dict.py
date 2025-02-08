# models/mutation_dict.py

from typing import Dict, Optional

class MutationDict:
    def __init__(self):
        self.entries: Dict[str, Dict] = {}

    def add_entry(self, entry_id: str, entry_data: Dict) -> None:
        """
        Add a new entry to the mutation dictionary.

        Args:
            entry_id (str): Unique identifier for the entry.
            entry_data (Dict): Data associated with the entry.
        """
        if entry_id in self.entries:
            raise ValueError(f"Entry with ID '{entry_id}' already exists.")
        self.entries[entry_id] = entry_data

    def get_entry(self, entry_id: str) -> Optional[Dict]:
        """
        Get an entry from the mutation dictionary by its ID.

        Args:
            entry_id (str): Unique identifier for the entry.

        Returns:
            Optional[Dict]: Entry data if found, None otherwise.
        """
        return self.entries.get(entry_id)

    def update_entry(self, entry_id: str, entry_data: Dict) -> None:
        """
        Update an existing entry in the mutation dictionary.

        Args:
            entry_id (str): Unique identifier for the entry.
            entry_data (Dict): Updated data for the entry.

        Raises:
            ValueError: If the entry with the given ID does not exist.
        """
        if entry_id not in self.entries:
            raise ValueError(f"Entry with ID '{entry_id}' does not exist.")
        self.entries[entry_id].update(entry_data)

    def remove_entry(self, entry_id: str) -> None:
        """
        Remove an entry from the mutation dictionary by its ID.

        Args:
            entry_id (str): Unique identifier for the entry.

        Raises:
            ValueError: If the entry with the given ID does not exist.
        """
        if entry_id not in self.entries:
            raise ValueError(f"Entry with ID '{entry_id}' does not exist.")
        del self.entries[entry_id]

    def __len__(self) -> int:
        """
        Get the number of entries in the mutation dictionary.

        Returns:
            int: Number of entries.
        """
        return len(self.entries)

    def __iter__(self):
        """
        Iterate over the entries in the mutation dictionary.

        Yields:
            Tuple[str, Dict]: Entry ID and associated data.
        """
        for entry_id, entry_data in self.entries.items():
            yield entry_id, entry_data

    def clear(self) -> None:
        """
        Clear all entries from the mutation dictionary.
        """
        self.entries.clear()