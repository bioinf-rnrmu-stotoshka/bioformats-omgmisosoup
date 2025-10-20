from abstract_reader import Reader
from abc import abstractmethod


class SequenceReader(Reader):
    """
    Abstract base class for biological sequence file readers.
    
    This class extends the Reader abstract class to provide common functionality
    for reading biological sequence files (FASTA, FASTQ, etc.) with sequence-specific
    methods and validation. Concrete subclasses must implement the abstract methods
    for format-specific reading and sequence retrieval.
    
    Attributes:
        file_path (str): Absolute path to the sequence file.
        
    Note:
        This class cannot be instantiated directly. Use format-specific subclasses
        like FastaReader or FastqReader instead.
    """


    def __init__(self, filename: str):
        """
        Initialize the sequence reader with a file path.
        
        Converts the provided filename to an absolute path and stores it
        for file reading operations.
        
        Args:
            filename (str): Path to the sequence file. Can be relative or absolute.
            
        Raises:
            ValueError: If filename is empty or None.
            
        Example:
            >>> reader = FastaReader("sequences.fasta")  # Concrete implementation
            >>> print(reader.file_path)
            '/absolute/path/to/sequences.fasta'
        """
        import os
        self.file_path = os.path.abspath(filename)
    

    @abstractmethod
    def read(self) -> dict:
        """
        Read and parse the sequence file.
        
        This method must be implemented by concrete subclasses to handle
        specific file formats. It should parse the file content and store
        sequences in an appropriate data structure.
        
        Returns:
            dict: Dictionary mapping sequence IDs to sequence strings.
            
        Raises:
            FileNotFoundError: If the specified file does not exist.
            ValueError: If the file format is invalid or contains malformed data.
            IOError: If there are issues reading the file.
            
        Note:
            Implementation must handle format-specific parsing logic and
            should call validate_sequence() to verify sequence content.
        """
        pass


    @abstractmethod
    def get_sequence(self, seq_id: str) -> str:
        """
        Retrieve a specific sequence by its ID.
        
        This method must be implemented by concrete subclasses to provide
        access to individual sequences stored in the reader.
        
        Args:
            seq_id (str): The unique identifier of the sequence.
            
        Returns:
            str: The biological sequence (nucleotide or protein) as a string.
            
        Raises:
            KeyError: If the sequence ID is not found in the stored data.
        """
        pass


    @abstractmethod
    def get_sequence_length(self, seq_id: str) -> int:
        """
        Get the length of a specific sequence.
        
        This method must be implemented by concrete subclasses to provide
        sequence length information.
        
        Args:
            seq_id (str): The unique identifier of the sequence.
            
        Returns:
            int: The length of the sequence in base pairs or amino acids.
            
        Raises:
            KeyError: If the sequence ID is not found in the stored data.
        """
        pass


    def validate_sequence(self, sequence: str) -> bool:
        """
        Validate nucleotide sequence format.
        
        Checks if the sequence contains only valid nucleotide characters
        (A, C, T, G, N in both uppercase and lowercase). N represents
        any nucleotide (unknown base).
        
        Args:
            sequence (str): Biological sequence string to validate.
            
        Returns:
            bool: 
                True - The sequence contains only valid nucleotide characters
                False - The sequence contains invalid characters
                
        Example:
            >>> reader.validate_sequence("ATCGNnatc")
            True
            >>> reader.validate_sequence("ATCGXnatc")  # 'X' is invalid
            False
            
        Note:
            This implementation validates DNA sequences. Subclasses may
            override this method for different sequence types (RNA, protein).
            The method is case-insensitive for nucleotide characters.
        """
        validate_set = set("ACTGactgNn")
        if set(sequence).issubset(validate_set):
            return True
        else:
            return False