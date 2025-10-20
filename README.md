# Bioinformatics File Readers

A Python library for reading and analyzing various bioinformatics file formats with comprehensive documentation and analysis capabilities.

## Supported Formats

- **FASTA** - Sequence files
- **FASTQ** - Sequence files with quality scores  
- **SAM** - Sequence Alignment/Map files
- **VCF** - Variant Call Format files

## Features

### Core Functionality
- **Abstract base classes** for consistent interface across all readers
- **Memory-efficient** processing using generators
- **Comprehensive validation** of file formats and sequences
- **Type annotations** for better code clarity and IDE support

### Analysis Capabilities
- **Quality score analysis** for FASTQ files
- **Variant statistics** for VCF files
- **Alignment metrics** for SAM files
- **Sequence statistics** for FASTA files

### Visualization
- **FastQC-style plots** for quality control
- **Statistical summaries** and reports
- **Genomic region filtering**

## Installation
```bash
git clone https://github.com/bioinf-rnrmu-stotoshka/bioformats-omgmisosoup.git
cd bioformats-omgmisosoup

python -m venv venv
source venv/bin/activate  # Linux/Mac
# или
venv\Scripts\activate     # Windows

# Установка в режиме разработки
pip install -e .
```

### Prerequisites
- Python 3.7+
- Required packages:

```bash
pip install pandas seaborn matplotlib numpy
```
### NOTE: For example of sam file DM solnyshko2009
  
