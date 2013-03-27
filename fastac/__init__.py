'''FastaC: A pure-python scripting-language and compiler for biological sequences.
by Cathal Garvey

FastaC stands for FASTA Compiler. It adopts and extends the FASTA file format,
commonly used in life sciences to represent biological sequences in a simple way,
adding comments, meta-comments, title-based JSON metadata, and scripting language
extensions enabling generation of complex sequences from precursors and procedurally
translated sequences at compile-time.

FastaC can be used simply to allow markup of normal FASTA files, in which case the
original (though nearly universally unsupported) ";" comment format is allowed
(and preserved on output as JSON metadata in the title, making it compatible with
other tools), as is an "ignored" comment using the "#" character, which is not
parsed or conserved by the compiler. This can allow normal Fasta files to be more
readable and informative to humans while compiling to FASTA that can be used easily
with the standard suite of applications in use already in bioinformatics.

Scripting extensions include:
* Construction of new blocks by concatenating or by performing transformations on
  previously-defined blocks, or on blocks imported from other files. Transformations
  may include reverse-complementation, translation, reverse translation, RNA<->DNA..
* Definition of "templates" which can be used as crude "functions" for encapsulating
  arbitrary sequences in a predefined context, for example an arbitrary CDS in a
  predefined expression frame consisting of a promoter, RBS, UTRs & terminator.
* Markup of blocks as "private", enabling internal definition of sub-blocks for
  scripting use only which are not themselves exported to compiled FASTA.
'''
__version__ = "0.1"
from fastac.compilefasta import FastaCompiler, FastaBlock, Macros, FastaError, FastaCompileError
