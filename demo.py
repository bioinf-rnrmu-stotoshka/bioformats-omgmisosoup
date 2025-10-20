"""
Simple Usage Examples for Bioinformatics Readers
"""

from classes.fasta_reader import FastaReader
from classes.fastq_reader import FastqReader
from classes.samreader import SAMReader
from classes.vcfreader import VCFProcessor

def simple_fasta_example():
    """Simple FASTA file reading example."""
    print("=== FASTA Example ===")
    reader = FastaReader(r"c:\Users\Asus\Desktop\miso\Sample.fasta")
    sequences = reader.read()
    
    print(f"Found {reader.count_sequences()} sequences")
    for seq_id, sequence in sequences.items():
        print(f"  {seq_id}: {len(sequence)} bp")
    
    return reader

def simple_fastq_example():
    """Simple FASTQ file reading example."""
    print("\n=== FASTQ Example ===")
    reader = FastqReader(r"c:\Users\Asus\Desktop\miso\ERR15416517.fastq")
    reader.read()
    
    print(f"Found {reader.count_sequences()} sequences")
    first_seq = list(reader.seq_dict.keys())[0]
    print(f"First sequence '{first_seq}':")
    print(f"  Length: {reader.get_sequence_length(first_seq)}")
    print(f"  Avg Quality: {reader.get_average_quality(first_seq)}")
    
    return reader

def simple_sam_example():
    """Simple SAM file reading example."""
    print("\n=== SAM Example ===")
    reader = SAMReader(r"c:\Users\Asus\Desktop\miso\exmpl.sam")
    
    print(f"Chromosomes: {reader.get_chromosomes()}")
    print(f"Total alignments: {reader.get_alignment_count()}")
    
    # Show first alignment
    for alignment in reader.read():
        print(f"First alignment: {alignment.qname} -> {alignment.rname}:{alignment.pos}")
        break
    
    reader.close()
    return reader

def simple_vcf_example():
    """Simple VCF file reading example."""
    print("\n=== VCF Example ===")
    processor = VCFProcessor(r"c:\Users\Asus\Desktop\miso\exmpl.vcf")
    processor.parse_file()
    
    print(f"Total variants: {processor.get_variant_count()}")
    
    # Show variant types
    stats = processor.get_variant_type_stats()
    print(f"Variant types: {stats}")
    
    processor.close()
    return processor

if __name__ == "__main__":
    # Run simple examples
    simple_fasta_example()
    simple_fastq_example() 
    simple_sam_example()
    simple_vcf_example()
    