from classes.sequence_reader import SequenceReader
from itertools import islice
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter


class FastqReader(SequenceReader):
    """
    A FASTQ file reader and analyzer for biological sequence data.
    
    This class extends SequenceReader to provide specialized functionality
    for reading and analyzing FASTQ format files, which contain DNA/RNA
    sequences along with quality scores.
    
    Attributes:
        total_sequences (int): Total number of sequences read from the file.
        total_length (int): Cumulative length of all sequences.
        seq_dict (dict): Dictionary mapping sequence IDs to sequence strings.
        seq_quality_dict (dict): Dictionary mapping sequence IDs to quality strings.
    
    Examples:
        >>> reader = FastqReader("sample.fastq")
        >>> reader.read()
        >>> print(f"Total sequences: {reader.count_sequences()}")
        >>> print(f"Average length: {reader.get_average_sequence_len()}")
    """

    def read(self):
        """
        Read and parse the FASTQ file, storing sequences and quality data.
        
        Processes the FASTQ file in 4-line chunks (header, sequence, separator, quality).
        Validates sequences and stores them along with their quality scores.
        Updates total sequence count and cumulative length statistics.
        
        Raises:
            FileNotFoundError: If the specified file path does not exist.
            ValueError: If the file format is invalid or sequences fail validation.
        """
        self.total_sequences = 0
        self.total_length = 0
        self.seq_dict = {}
        self.seq_quality_dict = {}

        with open(self.file_path) as file:
            while True:
                header = file.readline().strip()
                if not header:
                    break
                
                sequence = file.readline().strip()
                separator = file.readline().strip()
                quality = file.readline().strip()

                if self.validate_sequence(sequence):
                    seq_id = header[1:]
                    self.seq_dict[seq_id] = sequence
                    self.seq_quality_dict[seq_id] = quality
                    self.total_sequences += 1
                    self.total_length += len(sequence)
    

    def get_sequence(self, seq_id: str) -> str:
        """
        Retrieve a specific sequence by its ID.
        
        Args:
            seq_id (str): The unique identifier of the sequence (without the '@' prefix).
            
        Returns:
            str: The nucleotide sequence corresponding to the given ID.
            
        Raises:
            KeyError: If the sequence ID is not found in the stored data.
        """
        return self.seq_dict[seq_id]
    

    def get_sequence_length(self, seq_id: str) -> int:
        """
        Get the length of a specific sequence.
        
        Args:
            seq_id (str): The unique identifier of the sequence.
            
        Returns:
            int: The length of the sequence in base pairs.
            
        Raises:
            KeyError: If the sequence ID is not found in the stored data.
        """
        return len(self.seq_dict[seq_id])
    

    def count_sequences(self) -> int:
        """
        Count the total number of sequences loaded.
        
        Returns:
            int: The total number of sequences stored in the reader.
        """
        return len(self.seq_dict)
    

    def get_average_sequence_len(self) -> float:
        """
        Calculate the average length of all sequences.
        
        Returns:
            float: The average sequence length rounded to 2 decimal places.
            
        Note:
            Returns 0.0 if no sequences are loaded.
        """
        if self.total_sequences == 0:
            return 0.0
        return round(self.total_length / self.total_sequences, 2)
    

    def get_quality_scores(self, seq_id: str) -> list:
        """
        Convert quality string to numerical Phred quality scores.
        
        Args:
            seq_id (str): The unique identifier of the sequence.
            
        Returns:
            list: A list of integer Phred quality scores (ASCII value - 33).
            
        Raises:
            KeyError: If the sequence ID is not found in the stored data.
        """
        return [ord(char.encode(encoding="ascii")) - 33 for char in self.seq_quality_dict[seq_id]]
    

    def get_average_quality(self, seq_id: str) -> float:
        """
        Calculate the average quality score for a specific sequence.
        
        Args:
            seq_id (str): The unique identifier of the sequence.
            
        Returns:
            float: The average Phred quality score rounded to 2 decimal places.
            
        Raises:
            KeyError: If the sequence ID is not found in the stored data.
        """
        quality = self.get_quality_scores(seq_id)
        if not quality:
            return 0.0
        
        qual_sum = sum(quality)
        return round(qual_sum / len(quality), 2)
    

    def per_base_sequence_quality(self, sample_size=5000):
        """
        Generate a per-base sequence quality plot in FastQC style.
        
        Creates a boxplot showing quality score distribution across all sequence
        positions, with color-coded quality zones (green=good, orange=medium, red=poor).
        
        Args:
            sample_size (int, optional): Maximum number of sequences to sample for
                                        large files. Defaults to 5000.
                                        
        Raises:
            ValueError: If quality data is not loaded (read() not called).
        """
        if not hasattr(self, 'seq_quality_dict') or not self.seq_quality_dict:
            raise ValueError("Quality data not loaded. Call read() first.")
    
        # Sampling for large files
        qualities = list(self.seq_quality_dict.values())
        if len(qualities) > sample_size:
            import random
            qualities = random.sample(qualities, sample_size)
    
        max_len = max(len(qual) for qual in qualities)
    
        # Calculate quality statistics per position
        quality_stats = []
        for pos in range(max_len):
            pos_qualities = []
            for qual_str in qualities:
                if pos < len(qual_str):
                    pos_qualities.append(ord(qual_str[pos]) - 33)
        
            if pos_qualities:
                quality_stats.append({
                    'Position': pos + 1,
                    'Median': np.median(pos_qualities),
                    'Q25': np.percentile(pos_qualities, 25),
                    'Q75': np.percentile(pos_qualities, 75),
                    'Q10': np.percentile(pos_qualities, 10),
                    'Q90': np.percentile(pos_qualities, 90)
                })
    
        if not quality_stats:
            print("No quality data available for plotting")
            return
    
        quality_df = pd.DataFrame(quality_stats)

        fig, ax = plt.subplots(figsize=(12, 6))
    
        positions = quality_df['Position']
    
        # Create boxplot data
        stats = []
        for i, row in quality_df.iterrows():
            stats.append({
                'med': row['Median'],
                'q1': row['Q25'],
                'q3': row['Q75'], 
                'whislo': row['Q10'],
                'whishi': row['Q90'],
            })
    
        # Draw boxplot
        ax.bxp(stats, positions=range(len(positions)), 
               showfliers=False, patch_artist=True,
               boxprops=dict(facecolor='yellow', alpha=0.7),
               medianprops=dict(color='red', linewidth=2))
    
        # Quality zones
        ax.axhspan(28, 45, alpha=0.1, color='green')
        ax.axhspan(20, 28, alpha=0.1, color='orange')
        ax.axhspan(0, 20, alpha=0.1, color='red')
    
        ax.set_xlabel('Position in read (bp)')
        ax.set_ylabel('Quality scores')
        ax.set_title('Per base sequence quality')
        ax.set_xticks(range(len(positions)))
        ax.set_xticklabels(positions)
        ax.grid(True, alpha=0.3)
    
        plt.tight_layout()
        plt.show()
    
        # Memory cleanup
        del qualities, quality_stats


    def per_base_sequence_content(self, sample_size=5000):
        """
        Generate a per-base sequence content plot in FastQC style.
        
        Shows the percentage of each nucleotide (A, T, C, G) at each position
        across all sequences, with reference lines at 25% and 50%.
        
        Args:
            sample_size (int, optional): Maximum number of sequences to sample for
                                        large files. Defaults to 5000.
                                        
        Raises:
            ValueError: If sequence data is not loaded (read() not called).
        """
        if not hasattr(self, 'seq_dict') or not self.seq_dict:
            raise ValueError("Sequence data not loaded. Call read() first.")
    
        # Sampling
        sequences = list(self.seq_dict.values())
        if len(sequences) > sample_size:
            import random
            sequences = random.sample(sequences, sample_size)
    
        max_len = max(len(seq) for seq in sequences)
    
        # Count nucleotides per position
        base_counts = {base: np.zeros(max_len, dtype=np.int32) for base in 'ATCG'}
        total_counts = np.zeros(max_len, dtype=np.int32)
    
        for sequence in sequences:
            for pos, nucleotide in enumerate(sequence):
                if pos < max_len and nucleotide in base_counts:
                    base_counts[nucleotide][pos] += 1
                    total_counts[pos] += 1
    
        # Create DataFrame
        content_data = []
        for pos in range(max_len):
            if total_counts[pos] > 0:
                row = {'Position': pos + 1}
                for base in 'ATCG':
                    row[base] = (base_counts[base][pos] / total_counts[pos]) * 100
                content_data.append(row)
    
        content_df = pd.DataFrame(content_data)
        
        # Plot in FastQC style
        plt.figure(figsize=(12, 6))
    
        # Colors matching FastQC
        colors = {'A': 'red', 'T': 'blue', 'G': 'green', 'C': 'black'}
    
        # Draw lines for each nucleotide
        for base in 'ATCG':
            if base in content_df.columns:
                plt.plot(content_df['Position'], content_df[base], 
                        color=colors[base], linewidth=1.5, label=base)
    
        # FastQC-style settings
        plt.title('Per Base Sequence Content', fontsize=14, fontweight='bold')
        plt.xlabel('Position in read (bp)')
        plt.ylabel('Percentage')
    
        # Legend in upper right corner
        plt.legend(loc='upper right', frameon=True)
    
        # Grid only on Y axis
        plt.grid(True, axis='y', alpha=0.3, linestyle='-')
    
        # Reference lines at 25% and 50%
        plt.axhline(y=25, color='gray', linestyle=':', alpha=0.5, linewidth=0.8)
        plt.axhline(y=50, color='gray', linestyle=':', alpha=0.5, linewidth=0.8)
    
        # Axis ranges like FastQC
        plt.ylim(0, 100)
        if not content_df.empty:
            plt.xlim(1, min(len(content_df), 150))  # Typically show first 150 positions
    
        # Remove top and right borders
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
    
        plt.tight_layout()
        plt.show()
    
        # Memory cleanup
        del sequences, base_counts, total_counts, content_data


    def sequence_length_distribution(self):
        """
        Generate a sequence length distribution plot in FastQC style.
        
        Shows the distribution of sequence lengths across all loaded sequences
        using a histogram-style line plot.
        
        Raises:
            ValueError: If sequence data is not loaded (read() not called).
        """
        if not hasattr(self, 'seq_dict') or not self.seq_dict:
            raise ValueError("Sequence data not loaded. Call read() first.")
    
        lengths = [len(seq) for seq in self.seq_dict.values()]

        fig, ax = plt.subplots(figsize=(12, 6))
    
        counts, bin_edges = np.histogram(lengths, bins=30, density=False)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
        ax.plot(bin_centers, counts, color='blue', linewidth=2)
    
        ax.set_xlabel('Sequence Length (bp)')
        ax.set_title('Sequence Length Distribution')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)
    
        plt.tight_layout()
        plt.show()