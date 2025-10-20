from abc import ABC, abstractmethod


class Reader(ABC):
    """
    Abstract base class for biological data file readers.
    
    This class defines the interface for all biological sequence file readers,
    ensuring consistent implementation across different file formats (FASTA, FASTQ, etc.).
    Concrete subclasses must implement all abstract methods.
    
    Attributes:
        filename (str): Path to the file containing biological sequence data.
        
    Examples:
        class ConcreteReader(Reader):
            def __init__(self, filename):
                self.filename = filename
                
            def read(self):
                # Implementation for specific file format
                pass
    """
    

    @abstractmethod
    def __init__(self, filename):
        """
        Initialize the reader with a biological data file.
        
        Args:
            filename (str): Path to the biological data file to be read.
                           Supported formats depend on the concrete implementation.
                           
        Note:
            Concrete subclasses must call this constructor and store the filename.
        """
        pass


    @abstractmethod
    def read(self):
        """
        Read and parse the biological data file.
        
        This method must be implemented by concrete subclasses to handle
        specific file formats (FASTA, FASTQ, etc.). The implementation should
        parse the file content and store the sequences in an appropriate data structure.
        
        Raises:
            FileNotFoundError: If the specified file does not exist.
            ValueError: If the file format is invalid or contains malformed data.
            IOError: If there are issues reading the file.
            
        Note:
            Implementation details vary by file format. Subclasses should document
            format-specific requirements and behavior.
        """
        pass