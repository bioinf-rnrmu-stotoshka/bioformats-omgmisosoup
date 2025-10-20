import pandas as pd
import re
from typing import List, Dict, Iterator, Generator, Callable, Any
from genomicdatareader import GenomicDataReader


class Alignment:
    """Simple alignment data class"""
    def __init__(self, fields: List[str]):
        self.fields = fields
        self.qname = fields[0] if len(fields) > 0 else ""
        self.flag = int(fields[1]) if len(fields) > 1 else 0
        self.rname = fields[2] if len(fields) > 2 else ""
        self.pos = int(fields[3]) if len(fields) > 3 else 0
        self.mapq = int(fields[4]) if len(fields) > 4 else 0
        self.cigar = fields[5] if len(fields) > 5 else ""
        self.rnext = fields[6] if len(fields) > 6 else ""
        self.pnext = int(fields[7]) if len(fields) > 7 else 0
        self.tlen = int(fields[8]) if len(fields) > 8 else 0
        self.seq = fields[9] if len(fields) > 9 else ""
        self.qual = fields[10] if len(fields) > 10 else ""

    def __repr__(self):
        return f"Alignment({self.qname} -> {self.rname}:{self.pos})"


class SamHeader:
    """Simple SAM header data class"""
    def __init__(self, header_groups: Dict[str, List[str]]):
        self.header_groups = header_groups


class SAMReader(GenomicDataReader):
    """
    Memory-efficient SAM file reader using generators.
    Supports both DNA-seq and RNA-seq alignments with proper CIGAR parsing.
    """
    
    def __init__(self, filename: str):
        """Initialize with SAM file path."""
        super().__init__(filename)
        self._chromosomes = None
        self._file_handle = None
    
    def _read_lines(self) -> Generator[str, None, None]:
        """Generator that yields lines from SAM file."""
        with open(self.filename, 'r') as f:
            for line in f:
                yield line.strip()
    
    def read(self) -> Iterator[Alignment]:
        """Read alignments from SAM file."""
        for line in self._read_lines():
            if not line.startswith('@') and line:  # Skip header and empty lines
                fields = line.split('\t')
                if len(fields) >= 11:
                    yield self._parse_line(line)
    
    def _parse_line(self, line: str) -> Alignment:
        """Parse a SAM alignment line into Alignment object."""
        fields = line.split('\t')
        return Alignment(fields)
    
    def close(self):
        """Close any open file handles."""
        if self._file_handle:
            self._file_handle.close()
            self._file_handle = None
    
    def get_chromosomes(self) -> List[str]:
        """Get list of chromosomes from SQ headers."""
        if self._chromosomes is None:
            self._chromosomes = []
            header_groups = self.get_header_and_groups()
            for sq_line in header_groups['SQ']:
                # Extract chromosome name from @SQ line: @SQ\tSN:chr1\tLN:249250621
                for part in sq_line.split('\t'):
                    if part.startswith('SN:'):
                        self._chromosomes.append(part[3:])
                        break
        return self._chromosomes
    
    def get_reference_genome(self) -> str:
        """Get reference genome information from HD header."""
        header_groups = self.get_header_and_groups()
        for hd_line in header_groups['HD']:
            # Extract reference from @HD line if available
            for part in hd_line.split('\t'):
                if part.startswith('AS:'):
                    return part[3:]
        return "unknown"
    
    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """Validate if chromosome and position are valid."""
        if chrom not in self.get_chromosomes():
            return False
        if pos < 1:
            return False
        return True
    
    def get_header_and_groups(self) -> Dict[str, List[str]]:
        """
        Get header and information by header groups (@PG, @RG, etc.).
        
        Returns:
            Dictionary with header groups as keys and lines as values
        """
        groups = {'HD': [], 'SQ': [], 'RG': [], 'PG': [], 'CO': [], 'ALL': []}
        
        for line in self._read_lines():
            if line.startswith('@'):
                groups['ALL'].append(line)
                if line.startswith('@HD'):
                    groups['HD'].append(line)
                elif line.startswith('@SQ'):
                    groups['SQ'].append(line)
                elif line.startswith('@RG'):
                    groups['RG'].append(line)
                elif line.startswith('@PG'):
                    groups['PG'].append(line)
                elif line.startswith('@CO'):
                    groups['CO'].append(line)
            else:
                # Stop at first alignment line (headers are always first)
                break
                
        return groups
    
    def get_header(self) -> SamHeader:
        """Get SAM header as SamHeader object."""
        header_groups = self.get_header_and_groups()
        return SamHeader(header_groups)
    
    def read_alignments(self) -> List[Alignment]:
        """
        Read all alignments into memory.
        
        Returns:
            List of Alignment objects
        """
        return list(self.read())
    
    def iterate_alignments(self) -> Generator[List[str], None, None]:
        """
        Generator that yields alignments one by one as raw fields.
        
        Yields:
            List of alignment fields
        """
        for line in self._read_lines():
            if not line.startswith('@') and line:  # Skip header and empty lines
                fields = line.split('\t')
                if len(fields) >= 11:
                    yield fields
    
    def get_alignment_count(self) -> int:
        """
        Get total number of alignments using generator.
        
        Returns:
            Number of alignments
        """
        count = 0
        for _ in self.iterate_alignments():
            count += 1
        return count
    
    def filter_alignments(self, flag: int) -> List[Alignment]:
        """
        Filter alignments by flag.
        
        Args:
            flag: SAM flag to filter by
            
        Returns:
            List of filtered Alignment objects
        """
        filtered = []
        for alignment in self.iterate_alignments():
            aln_flag = int(alignment[1])
            if aln_flag & flag:
                filtered.append(Alignment(alignment))
        return filtered
    
    def calculate_coverage(self, chrom: str) -> Dict[int, int]:
        """
        Calculate coverage for a chromosome.
        
        Args:
            chrom: Chromosome name
            
        Returns:
            Dictionary with position -> coverage
        """
        coverage = {}
        for alignment in self.iterate_alignments():
            if alignment[2] == chrom:  # RNAME field
                try:
                    pos = int(alignment[3])  # POS field
                    cigar = alignment[5]     # CIGAR field
                    
                    # Simple coverage calculation - extend coverage by read length
                    read_length = len(alignment[9]) if alignment[9] != '*' else 0
                    for i in range(pos, pos + read_length):
                        coverage[i] = coverage.get(i, 0) + 1
                except (ValueError, IndexError):
                    continue
        return coverage
    
    def get_chromosome_stats(self) -> pd.DataFrame:
        """
        Get alignment counts by chromosome using pandas and generators.
        
        Returns:
            DataFrame with chromosome counts
        """
        from collections import defaultdict
        chrom_counts = defaultdict(int)
        
        for alignment in self.iterate_alignments():
            chrom = alignment[2]  # RNAME field
            if chrom != '*':  # Skip unmapped
                chrom_counts[chrom] += 1
        
        df = pd.DataFrame(list(chrom_counts.items()), 
                         columns=['Chromosome', 'Count'])
        return df.sort_values('Count', ascending=False)
    
    def _parse_cigar(self, cigar: str) -> List[tuple]:
        """
        Parse CIGAR string into list of (length, operation) tuples.
        
        Args:
            cigar: CIGAR string (e.g., "43S13M6872N45M")
            
        Returns:
            List of tuples [(length, operation), ...]
        """
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
    
    def _alignment_overlaps_region(self, pos: int, cigar: str, start: int, end: int) -> bool:
        """
        Check if alignment overlaps with genomic region considering CIGAR operations.
        
        Args:
            pos: Start position of alignment
            cigar: CIGAR string
            start: Region start position
            end: Region end position
            
        Returns:
            True if alignment overlaps the region, False otherwise
        """
        current_pos = pos
        cigar_ops = self._parse_cigar(cigar)
        
        for length, op in cigar_ops:
            if op in ['M', '=', 'X', 'D', 'N']:  # Operations that consume reference
                align_start = current_pos
                align_end = current_pos + length - 1
                
                # Check if this segment overlaps with the target region
                if max(align_start, start) <= min(align_end, end):
                    return True
                
                current_pos += length
            elif op in ['S', 'H', 'I']:  # Operations that don't consume reference
                continue
        
        return False
    
    def iterate_alignments_in_region(self, chrom: str, start: int, end: int) -> Generator[List[str], None, None]:
        """
        Generator for alignments in a genomic region with CIGAR-aware filtering.
        
        Args:
            chrom: Chromosome name
            start: Start position (1-based)
            end: End position (1-based)
            
        Yields:
            Alignments that overlap the specified region
        """
        if not self.validate_coordinate(chrom, start):
            raise ValueError(f"Invalid coordinate: {chrom}:{start}")
            
        for alignment in self.iterate_alignments():
            if alignment[2] == chrom:  # Check chromosome
                try:
                    pos = int(alignment[3])  # POS field
                    cigar = alignment[5]     # CIGAR field
                    
                    # Check if alignment overlaps region considering CIGAR
                    if self._alignment_overlaps_region(pos, cigar, start, end):
                        yield alignment
                except (ValueError, IndexError):
                    continue  # Skip if position is not a number or CIGAR missing
    
    def get_alignments_in_region(self, chrom: str, start: int, end: int) -> List[List[str]]:
        """
        Get all alignments in a genomic region (convenience method).
        
        Args:
            chrom: Chromosome name
            start: Start position
            end: End position
            
        Returns:
            List of alignments in the region
        """
        return list(self.iterate_alignments_in_region(chrom, start, end))


# Example usage with generators
if __name__ == "__main__":
    sam_reader = SAMReader(r"c:\Users\Asus\Desktop\sasam\eeebat.sam")
    
    # Test abstract class methods
    print("Filename:", sam_reader.filename)
    print("Chromosomes:", sam_reader.get_chromosomes())
    print("Reference genome:", sam_reader.get_reference_genome())
    print("Coordinate validation chr1:1000:", sam_reader.validate_coordinate("chr1", 1000))
    
    # Test reading using abstract interface
    print("\nReading alignments using abstract read() method:")
    alignment_count = 0
    for alignment in sam_reader.read():
        alignment_count += 1
        if alignment_count <= 3:  # Show first 3
            print(f"  {alignment}")
    
    print(f"Total alignments read: {alignment_count}")
    
    # Close reader
    sam_reader.close()