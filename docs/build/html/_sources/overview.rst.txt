Overview
========

This package provides a set of classes for reading and analyzing bioinformatics files:

- **FASTA**: Reading sequence files
- **FASTQ**: Reading sequence files with quality scores
- **SAM**: Reading alignment files
- **VCF**: Reading genomic variant files

Architecture
------------

The following diagram shows the class hierarchy and relationships:

.. figure:: _static/uml.png
   :width: 100%
   :alt: Bioinformatics Readers UML Class Diagram
   :align: center
   :class: with-shadow

   Figure 1: UML Class Diagram of Bioinformatics Readers


The package follows an object-oriented design with abstract base classes:

.. inheritance-diagram:: abstract_reader.Reader sequence_reader.SequenceReader genomicdatareader.GenomicDataReader
   :parts: 1

Main Components
---------------

Abstract Base Classes
^^^^^^^^^^^^^^^^^^^^^

* :class:`abstract_reader.Reader` - Base reader interface
* :class:`sequence_reader.SequenceReader` - Base sequence reader with validation
* :class:`genomicdatareader.GenomicDataReader` - Base genomic data reader

Concrete Implementations
^^^^^^^^^^^^^^^^^^^^^^^^

* :class:`fasta_reader.FastaReader` - FASTA file reader
* :class:`fastq_reader.FastqReader` - FASTQ file reader with quality analysis
* :class:`samreader.SAMReader` - SAM alignment file reader
* :class:`vcfreader.VCFProcessor` - VCF variant file processor

Usage Examples
--------------

FASTA Files
^^^^^^^^^^^

.. code-block:: python

   from fasta_reader import FastaReader
   
   reader = FastaReader("sequences.fasta")
   reader.read()
   sequence = reader.get_sequence("seq1")
   print(f"Sequence length: {reader.get_sequence_length('seq1')}")

FASTQ Files
^^^^^^^^^^^

.. code-block:: python

   from fastq_reader import FastqReader
   
   reader = FastqReader("sample.fastq")
   reader.read()
   reader.per_base_sequence_quality()
   print(f"Average quality: {reader.get_average_quality('read1')}")

SAM Files
^^^^^^^^^

.. code-block:: python

   from samreader import SAMReader
   
   reader = SAMReader("alignment.sam")
   chromosomes = reader.get_chromosomes()
   for alignment in reader.read():
       print(alignment)

VCF Files
^^^^^^^^^

.. code-block:: python

   from vcfreader import VCFProcessor
   
   processor = VCFProcessor("variants.vcf")
   processor.parse_file()
   report = processor.generate_summary_report()
   high_quality_variants = processor.filter_by_quality(20)