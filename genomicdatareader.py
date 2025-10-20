from abc import ABC, abstractmethod
from typing import List, Dict, Any, Iterator
from abstract_reader import Reader


class GenomicDataReader(Reader):
    """
    Abstract base class for genomic data readers
    Inherits from Reader and adds genomic-specific functionality
    """
    
    def __init__(self, filename: str):
        """Initialize with genomic data file."""
        super().__init__(filename)
    
    @abstractmethod
    def get_chromosomes(self) -> List[str]:
        """Get list of chromosomes in the data"""
        pass
    
    @abstractmethod
    def get_reference_genome(self) -> str:
        """Get reference genome information"""
        pass
    
    @abstractmethod
    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """Validate if chromosome and position are valid"""
        pass
    
    @abstractmethod
    def read(self) -> Iterator[Any]:
        """Read genomic records from file."""
        pass
    
    @abstractmethod
    def _parse_line(self, line: str) -> Any:
        """Parse a single line into a genomic Record object."""
        pass