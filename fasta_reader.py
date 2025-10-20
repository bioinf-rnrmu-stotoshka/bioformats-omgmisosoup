from sequence_reader import SequenceReader


class FastaReader(SequenceReader):
    """
    A FASTA file reader and analyzer for biological sequence data.
    
    This class extends SequenceReader to provide specialized functionality
    for reading and analyzing FASTA format files, which contain DNA/RNA
    sequences with header lines starting with '>'.
    
    Attributes:
        seq_dict (dict): Dictionary mapping sequence IDs to sequence strings.
        file_path (str): Path to the FASTA file being read.
    
    Examples:
        >>> reader = FastaReader("sequences.fasta")
        >>> reader.read()
        >>> print(f"Total sequences: {reader.count_sequences()}")
        >>> print(f"Average length: {reader.get_average_sequence_len()}")
        >>> sequence = reader.get_sequence("seq1")
    """

    def read(self) -> dict:
        """
        Read and parse the FASTA file, storing sequences in a dictionary.
        
        Processes the FASTA file line by line, handling multi-line sequences
        and validating both sequence headers and content. Sequences are stored
        in a dictionary with the header (without '>') as the key.
        
        Returns:
            dict: Dictionary containing sequence IDs as keys and sequences as values.
            
        Raises:
            ValueError: If file path is not provided and user input is empty.
            ValueError: If an empty sequence tag (header) is encountered.
            ValueError: If an empty sequence is found.
            ValueError: If sequence validation fails (invalid characters).
            FileNotFoundError: If the specified file path does not exist.
            IOError: If there are issues reading the file.
            
        Note:
            If file_path is not set, will prompt user for input.
            Sequences are validated using the parent class's validate_sequence method.
            Multi-line sequences are concatenated into single sequence strings.
        """
        self.seq_dict = {}
        if self.file_path is None:
            self.file_path = input("Enter the path to the FASTA-file: ")
        seq = ""
        tag = None
        with open(self.file_path) as file:
            for line in file:
                if not line:
                    continue
                line = line.strip()
                if line.startswith(">"):
                    if tag is not None:
                        if not seq:
                            raise ValueError("Empty sequence is found")
                        if self.validate_sequence(seq):
                            self.seq_dict[tag] = seq
                        else:
                            raise ValueError("Incorrect sequence format")
                        seq = ""
                    tag = line[1:]
                    if not tag:
                        raise ValueError("Empty tag is found")
                else:
                    seq += line
            if tag is not None:
                if not seq:
                    raise ValueError("Empty sequence is found")
                if self.validate_sequence(seq):
                    self.seq_dict[tag] = seq
                else:
                    raise ValueError("Incorrect sequence found")
        
        return self.seq_dict
    
    def get_sequence(self, seq_id: str) -> str:
        """
        Retrieve a specific sequence by its ID.
        
        Args:
            seq_id (str): The unique identifier of the sequence (header without '>').
            
        Returns:
            str: The nucleotide/protein sequence corresponding to the given ID.
            
        Raises:
            KeyError: If the sequence ID is not found in the stored data.
            
        Example:
            >>> reader.get_sequence("NP_001234")
            'ATGCTAGCTAGCTACGATCGATCGATCG...'
        """
        return self.seq_dict[seq_id]
    

    def get_sequence_length(self, seq_id: str) -> int:
        """
        Get the length of a specific sequence.
        
        Args:
            seq_id (str): The unique identifier of the sequence.
            
        Returns:
            int: The length of the sequence in base pairs or amino acids.
            
        Raises:
            KeyError: If the sequence ID is not found in the stored data.
            
        Example:
            >>> reader.get_sequence_length("NP_001234")
            150
        """
        return len(self.seq_dict[seq_id])
    

    def count_sequences(self) -> int:
        """
        Count the total number of sequences loaded.
        
        Returns:
            int: The total number of sequences stored in the reader.
            
        Example:
            >>> reader.count_sequences()
            25
        """
        return len(self.seq_dict)
    

    def get_average_sequence_len(self) -> float:
        """
        Calculate the average length of all sequences.
        
        Returns:
            float: The average sequence length rounded to 2 decimal places.
            
        Note:
            Returns 0.0 if no sequences are loaded.
            
        Example:
            >>> reader.get_average_sequence_len()
            245.67
        """
        if not self.seq_dict:
            return 0.0
            
        len_sum = 0
        for seq in self.seq_dict:
            len_sum += self.get_sequence_length(seq)

        return round(len_sum / len(self.seq_dict), 2)