import pandas as pd
import re
from typing import List, Dict, Any, Iterator
from genomicdatareader import GenomicDataReader


class Variant:
    """Simple variant data class"""
    def __init__(self, chrom: str, pos: int, ref: str, alt: str, qual: float = None, 
                 variant_id: str = ".", filter_status: str = ".", info: Dict = None):
        self.chrom = chrom
        self.pos = pos
        self.id = variant_id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter_status
        self.info = info or {}

    def __repr__(self):
        return f"Variant({self.chrom}:{self.pos} {self.ref}->{self.alt})"


class Genotype:
    """Simple genotype data class"""
    def __init__(self, sample: str, alleles: List[str], phased: bool = False):
        self.sample = sample
        self.alleles = alleles
        self.phased = phased

    def __repr__(self):
        phase_char = "|" if self.phased else "/"
        return f"Genotype({self.sample}: {phase_char.join(self.alleles)})"


class VcfHeader:
    """Simple VCF header data class"""
    def __init__(self, metadata: Dict, samples: List[str]):
        self.metadata = metadata
        self.samples = samples


class VCFProcessor(GenomicDataReader):
    """
    Unified class for working with VCF files
    Combines functionality for reading, analyzing and generating reports
    """

    def __init__(self, vcf_file):
        """
        Initialize the class
        
        Args:
            vcf_file (str): Path to VCF file
        """
        super().__init__(vcf_file)
        self.metadata = {}
        self.header = []
        self.samples = []
        self.variants_df = None
        self._chromosomes = None
        self._file_handle = None

    def read(self) -> Iterator[Variant]:
        """Read variants from VCF file."""
        if self.variants_df is None:
            self.parse_file()
            
        for _, row in self.variants_df.iterrows():
            yield self._parse_line_from_row(row)
    
    def _parse_line(self, line: str) -> Variant:
        """Parse a VCF variant line into Variant object."""
        fields = line.split('\t')
        if len(fields) < 8:
            raise ValueError(f"Invalid VCF line: {line}")
            
        return Variant(
            chrom=fields[0],
            pos=int(fields[1]),
            ref=fields[3],
            alt=fields[4],
            qual=self._parse_qual(fields[5]),
            variant_id=fields[2],
            filter_status=fields[6],
            info=self._parse_info_field(fields[7])
        )
    
    def _parse_line_from_row(self, row) -> Variant:
        """Parse a row from DataFrame into Variant object."""
        info_dict = {k: v for k, v in row.items() if k not in ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter']}
        return Variant(
            chrom=row['chrom'],
            pos=row['pos'],
            ref=row['ref'],
            alt=row['alt'],
            qual=row.get('qual'),
            variant_id=row.get('id', '.'),
            filter_status=row.get('filter', '.'),
            info=info_dict
        )
    
    def close(self):
        """Close any open file handles."""
        if self._file_handle:
            self._file_handle.close()
            self._file_handle = None

    def get_chromosomes(self) -> List[str]:
        """Get list of chromosomes from variants data."""
        if self._chromosomes is None:
            if self.variants_df is not None and not self.variants_df.empty:
                self._chromosomes = self.variants_df['chrom'].unique().tolist()
            else:
                self._chromosomes = []
        return self._chromosomes
    
    def get_reference_genome(self) -> str:
        """Get reference genome information from metadata."""
        return self.metadata.get('fileformat', 'unknown')
    
    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """Validate if chromosome and position are valid."""
        if chrom not in self.get_chromosomes():
            return False
        if pos < 1:
            return False
        return True

    def _parse_metadata(self, metadata_lines):
        """
        Parse VCF file metadata
        
        Args:
            metadata_lines (list): List of metadata lines
        """
        self.metadata = {'fileformat': '', 'info': {}, 'format': {}, 'filter': {}, 'other': {}}

        for line in metadata_lines:
            if line.startswith('##fileformat='):
                self.metadata['fileformat'] = line.split('=', 1)[1]
            elif line.startswith('##INFO='):
                # Parse INFO field: ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
                match = re.match(r'##INFO=<ID=([^,]+),.*,Description="([^"]+)">', line)
                if match:
                    self.metadata['info'][match.group(1)] = match.group(2)
            elif line.startswith('##FORMAT='):
                # Parse FORMAT field: ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                match = re.match(r'##FORMAT=<ID=([^,]+),.*,Description="([^"]+)">', line)
                if match:
                    self.metadata['format'][match.group(1)] = match.group(2)
            elif line.startswith('##FILTER='):
                # Parse FILTER field: ##FILTER=<ID=PASS,Description="All filters passed">
                match = re.match(r'##FILTER=<ID=([^,]+),Description="([^"]+)">', line)
                if match:
                    self.metadata['filter'][match.group(1)] = match.group(2)
            elif '=' in line[2:]:
                key, value = line[2:].split('=', 1)
                self.metadata['other'][key] = value

    def _parse_info_field(self, info_str):
        """
        Parse INFO field into dictionary
        
        Args:
            info_str (str): INFO field string
            
        Returns:
            dict: Dictionary with parsed values
        """
        info_dict = {}
        if info_str != '.':
            for item in info_str.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    try:
                        if '.' in value:
                            value = float(value)
                        else:
                            value = int(value)
                    except:
                        pass
                    info_dict[key.lower()] = value
                else:
                    info_dict[item.lower()] = True
        return info_dict

    def _parse_qual(self, qual_str):
        """
        Parse quality field
        
        Args:
            qual_str (str): Quality string
            
        Returns:
            float or None: Numeric quality value or None
        """
        if qual_str != '.':
            try:
                return float(qual_str)
            except:
                return None
        return None

    def parse_file(self):
        """
        Read and parse VCF file
        Creates pandas dataframe with variants
        """
        print(f"Reading VCF file: {self.filename}")

        metadata_lines = []
        header_line = ""
        variant_lines = []

        with open(self.filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('##'):
                    metadata_lines.append(line)
                elif line.startswith('#CHROM'):
                    header_line = line
                else:
                    variant_lines.append(line)

        self._parse_metadata(metadata_lines)
        self._parse_header(header_line)

        if variant_lines:
            self.variants_df = self._create_dataframe(variant_lines)
            print(f"Total variants: {len(self.variants_df)}")
        else:
            print("No variants found")

    def _parse_header(self, header_line):
        """
        Parse VCF header line
        
        Args:
            header_line (str): Header line
        """
        if header_line:
           self.header = header_line[1:].split('\t')
           if len(self.header) > 9:
                self.samples = self.header[9:]
                print(f"Number of samples: {len(self.samples)}")

    def _create_dataframe(self, variant_lines):
        """
        Create pandas dataframe from variant lines
        
        Args:
            variant_lines (list): List of variant lines
            
        Returns:
            pandas.DataFrame: Dataframe with variants
        """
        data = []
        for line in variant_lines:
            fields = line.split('\t')
            if len(fields) >= 8:
                variant_data = {
                    'chrom': fields[0],
                    'pos': int(fields[1]),
                    'id': fields[2],
                    'ref': fields[3],
                    'alt': fields[4],
                    'qual': self._parse_qual(fields[5]),
                    'filter': fields[6],
                    'info': fields[7]
                }

                # Add information from INFO field
                info_dict = self._parse_info_field(fields[7])
                variant_data.update(info_dict)

                data.append(variant_data)

        return pd.DataFrame(data)

    def get_header(self) -> VcfHeader:
        """Get VCF header as VcfHeader object."""
        return VcfHeader(self.metadata, self.samples)

    def read_variants(self) -> List[Variant]:
        """
        Read all variants into memory as Variant objects.
        
        Returns:
            List of Variant objects
        """
        return list(self.read())

    def filter_by_quality(self, min_qual: float) -> List[Variant]:
        """
        Filter variants by quality threshold.
        
        Args:
            min_qual: Minimum quality score
            
        Returns:
            List of filtered Variant objects
        """
        all_variants = self.read_variants()
        return [v for v in all_variants if v.qual is not None and v.qual >= min_qual]

    def get_genotype(self, sample: str, variant: Variant) -> Genotype:
        """
        Get genotype for a specific sample and variant.
        Note: This is a simplified implementation.
        
        Args:
            sample: Sample name
            variant: Variant object
            
        Returns:
            Genotype object
        """
        # This would need to be implemented based on your specific VCF format
        # For now, returning a placeholder implementation
        return Genotype(sample=sample, alleles=[variant.ref, variant.alt], phased=False)


    def get_header_info(self, group=None):
        """
        Get header and metadata information
        
        Args:
            group (str, optional): Specific metadata group ('info', 'format', 'filter')
            
        Returns:
            dict: Dictionary with information about requested group
        """
        if not self.metadata:
            print("No metadata available")
            return {}

        if group:
            return self.metadata.get(group.lower(), {})
        else:
            return {
                'fileformat': self.metadata['fileformat'],
                'samples': self.samples,
                'header_columns': self.header
            }

    def get_variant_count(self):
        """
        Get total number of variants
        
        Returns:
            int: Number of variants in the file
        """
        if self.variants_df is None:
            print("Data not loaded")
            return 0

        count = len(self.variants_df)
        print(f"Variant count: {count}")
        return count

    def get_alignment_statistics(self, region=None):
        """
        Get alignment statistics by regions
        
        Args:
            region (str, optional): Specific chromosome for filtering
            
        Returns:
            pandas.DataFrame: Dataframe with statistics
        """
        if self.variants_df is None:
            print("Data not loaded")
            return pd.DataFrame()

        # Polymorphism - different behavior depending on region
        if region:
            data = self.variants_df[self.variants_df['chrom'] == region]
            if data.empty:
                print(f"No variants found for region {region}")
                return pd.DataFrame()
        else:
            data = self.variants_df

        # Group by chromosome and calculate statistics
        stats = data.groupby('chrom').agg({
            'pos': 'count',  # variant count
            'dp': ['mean', 'median'] if 'dp' in data.columns else 'count'
        }).round(2)

        if 'dp' in data.columns:
            stats.columns = ['variant_count', 'mean_dp', 'median_dp']
        else:
            stats.columns = ['variant_count']

        return stats.reset_index()

    def get_variants_in_region(self, chrom, start, end):
        """
        Find variants in specified genomic region
        
        Args:
            chrom (str): Chromosome
            start (int): Start position
            end (int): End position
            
        Returns:
            pandas.DataFrame: Dataframe with variants in the region
        """
        if self.variants_df is None:
            print("Data not loaded")
            return pd.DataFrame()

        # Filter variants by position
        variants_in_region = self.variants_df[
            (self.variants_df['chrom'] == chrom) &
            (self.variants_df['pos'] >= start) &
            (self.variants_df['pos'] <= end)
        ]

        print(f"Found variants in region {chrom}:{start}-{end}: {len(variants_in_region)}")
        return variants_in_region

    def get_variant_type_stats(self):
        """
        Get statistics by variant types
        
        Returns:
            dict: Dictionary with variant counts by type
        """
        if self.variants_df is None:
            print("No data available")
            return {}

        stats = {'snp': 0, 'indel': 0, 'other': 0}

        for _, variant in self.variants_df.iterrows():
            ref_len = len(variant['ref'])
            alt_len = len(variant['alt'])

            if ref_len == 1 and alt_len == 1:
                stats['snp'] += 1
            elif ref_len != alt_len:
                stats['indel'] += 1
            else:
                stats['other'] += 1

        return stats

    def generate_summary_report(self):
        """
        Generate summary report for VCF file
        
        Returns:
            dict: Dictionary with overall statistics
        """
        if self.variants_df is None:
            print("No data available")
            return {}

        report = {
            'total_variants': len(self.variants_df),
            'chromosomes': self.variants_df['chrom'].unique().tolist(),
            'variant_types': self.get_variant_type_stats(),
            'region_statistics': self.get_alignment_statistics().to_dict('records')
        }

        # Add quality statistics if available
        if 'qual' in self.variants_df.columns and self.variants_df['qual'].notna().any():
            report['quality'] = {
                'mean': round(self.variants_df['qual'].mean(), 2),
                'median': round(self.variants_df['qual'].median(), 2)
            }

        # Add depth statistics if available
        if 'dp' in self.variants_df.columns and self.variants_df['dp'].notna().any():
            report['coverage_depth'] = {
                'mean': round(self.variants_df['dp'].mean(), 2),
                'median': round(self.variants_df['dp'].median(), 2)
            }

        return report

    def save_variants_to_file(self, variants_df, output_file):
        """
        Save variants to file
        
        Args:
            variants_df (pandas.DataFrame): Dataframe with variants
            output_file (str): Output file path
        """
        if variants_df.empty:
            print("No data to save")
            return

        variants_df.to_csv(output_file, sep='\t', index=False)


def main():
    vcf_file = "/content/S1.haplotypecaller.filtered.phased.csq.vcf"
    processor = VCFProcessor(vcf_file)
    processor.parse_file()

    if processor.variants_df is None:
        print("No data available")
        return

    # Test abstract class methods
    print("Chromosomes:", processor.get_chromosomes())
    print("Reference genome:", processor.get_reference_genome())
    print("Coordinate validation chr1:69270:", processor.validate_coordinate("chr1", 69270))

    # Header information
    print("Header Information:")
    header_info = processor.get_header_info()
    print(f"Format: {header_info.get('fileformat', '')}")
    print(f"Samples: {header_info.get('samples', [])}")

    # Variant count
    print("Variant Count:")
    count = processor.get_variant_count()

    # Test reading variants as objects
    print("Reading variants as objects:")
    variants = processor.read_variants()
    print(f"Read {len(variants)} variant objects")

    # Test quality filtering
    print("Quality filtering (QUAL >= 20):")
    high_qual_variants = processor.filter_by_quality(20)
    print(f"Found {len(high_qual_variants)} high quality variants")

    # Alignment statistics
    print("Alignment Statistics:")
    stats = processor.get_alignment_statistics()
    if not stats.empty:
        print(stats.to_string(index=False))

    # Region search
    print("Searching variants in region:")
    region_variants = processor.get_variants_in_region("chr1", 69270, 200000)
    if not region_variants.empty:
        print("First 5 variants:")
        cols = ['chrom', 'pos', 'ref', 'alt', 'qual']
        available_cols = [c for c in cols if c in region_variants.columns]
        print(region_variants[available_cols].head().to_string(index=False))

    report = processor.generate_summary_report()
    print(f"Chromosomes: {report.get('chromosomes', [])}")
    print(f"Variant types: {report.get('variant_types', {})}")


if __name__ == "__main__":
    main()