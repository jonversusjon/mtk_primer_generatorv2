# class MutationDict(dict):
#     """
#     A dictionary class for storing mutation entries that inherits from dict.
#     This gives us all standard dictionary functionality plus our custom methods.
#     """
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)

#     def add_entry(self, entry_id: str, entry_data: dict) -> None:
#         """
#         Add a new entry to the mutation dictionary.

#         Args:
#             entry_id (str): Unique identifier for the entry.
#             entry_data (dict): Data associated with the entry.
#         """
#         if entry_id in self:
#             raise ValueError(f"Entry with ID '{entry_id}' already exists.")
#         self[entry_id] = entry_data

#     def get_entry(self, entry_id: str) -> dict:
#         """
#         Get an entry from the mutation dictionary by its ID.

#         Args:
#             entry_id (str): Unique identifier for the entry.

#         Returns:
#             dict: Entry data if found, None otherwise.
#         """
#         return self.get(entry_id)

#     def update_entry(self, entry_id: str, entry_data: dict) -> None:
#         """
#         Update an existing entry in the mutation dictionary.

#         Args:
#             entry_id (str): Unique identifier for the entry.
#             entry_data (dict): Updated data for the entry.

#         Raises:
#             ValueError: If the entry with the given ID does not exist.
#         """
#         if entry_id not in self:
#             raise ValueError(f"Entry with ID '{entry_id}' does not exist.")
#         self[entry_id].update(entry_data)

#     def remove_entry(self, entry_id: str) -> None:
#         """
#         Remove an entry from the mutation dictionary by its ID.

#         Args:
#             entry_id (str): Unique identifier for the entry.

#         Raises:
#             ValueError: If the entry with the given ID does not exist.
#         """
#         if entry_id not in self:
#             raise ValueError(f"Entry with ID '{entry_id}' does not exist.")
#         del self[entry_id]

#     # No need to define __len__ or __iter__ as these come from dict
    
#     def clear(self) -> None:
#         """
#         Clear all entries from the mutation dictionary.
#         """
#         super().clear()