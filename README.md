# FASTAC - A FASTA Compiler
(C) Cathal Garvey, Licensed under the GNU Affero General Public Licence v3

Twitter: @onetruecathal

Neglected Blog: http://www.indiebiotech.com

## What is this?
This is a "compiler" for FASTA files which supports markup in the original FASTA style, using semicolons, and additionally supports a number of bash-scripting-style extensions including "macros" such as sequence inclusion from "libraries", translation, basic reverse translation, and reverse complementation. Additionally, it supports metadata definition whether in the title as a large block (allowing metadata to be used in a manner that does not interfere with other FASTA-compatible tools), or inline using semicolon markup. Finally, 'throwaway' comments are permitted using the "#" character, which are ignored during compilation and are not entered as metadata.

In other words, it's "yet another synthetic biology compiler".

## Huh?
The easiest way to grok my jive is to open up a bash prompt, move into the tests directory, and call "../fastac.py testcase.fasta", then compare the raw file to the output you see on the terminal. The raw files, including the two demo "libraries", make heavy use of the scripting extensions to define compound or transformed blocks, which are then imported and further transformed in the "testcase.fasta" file. The output, however, is standard FASTA, formed to 50 characters per line with any remaining metadata compressed into the title line.

## How do I use this?
You'll need to acquaint yourself with a terminal prompt such as Bash/Sh or (if you use Windows) cmd.exe. I leave that to you to figure out.

Then, just call this using python3 (Linux users can just call the script directly after using chmod to set the executable flag):
$ python3 fastac.py [fastafile]

For additional options, mostly at this point relating to the style of the output, try:
$ python3 fastac.py --help

Right now, output contains metadata in the title by default, such as the (inferred or user-defined) sequence type and any inline ";" comments: I'll add an option to disable this in future.

## How do I extend this?
To write your own "scripts", simply use valid multi-fasta files as "libraries", and then "include" sequences from these libraries using similar commands to those you see in the testcase.fasta file.

"Macro" commands are invoked on their own lines within a fasta block with a "$" character before the macro. Macro commands have a bash-style calling syntax, with positional and optional arguments. Examples of usage are given in the testcase.fasta file.

When a reference to a Fasta block is necessary (as with include, for example), use the title minus any metadata JSON blocks. So, a block with a title consisting of "> LacI Amino Acid Sequence {"type":"aminos"}" should be referred to as "LacI Amino Acid Sequence", with quotes. If the title has no space characters (recommended), then quotes are not needed.

Available macros as of first upload are:
include [--lib libraryfile] fasta_block_by_title
complement [--lib libraryfile] fasta_block_by_title
translate [--lib libraryfile] [--table table_name] fasta_block_by_title
dumb_backtranslate [--lib libraryfile] [--table table_name] fasta_block_by_title

Adding new macros is ~easy: Just define a function that takes a list of arguments as returned by the shlex.split() function in the Python standard library, and a FastaCompiler object (which defines the current scope for the macro, allowing libraries to recurse).

Your function should return the results of whatever transforms it has been called to perform for direct inclusion in Fasta blocks in which it is called. Then add your new function to the Macros dictionary with the name it will be called by, preferably the same name as the function itself.

The easiest way to handle lists of arguments as returned by shlex.split is to define an argparse.ArgumentParser instance and call "add_argument()" on the argparse instance to add each argument you expect or support; check the Python Standard Library for help on how to use argparse, or mimic what I've done with the builtin Macros. Then tell your parser to parse the args list, and use the returned namespace object to get the passed argument values.

#### That's horrible
Get over it.

## Why?
There are plenty of projects to make DNA "compilers" that operate on a GUI-only basis, emphasising drag-and-drop design of DNA "parts". However, most of these are closed source, and most also function essentially as a pretty veneer on copy/paste operations. I've endeavoured to create in fastac a simple but extensible scripting language which can be enhanced easily with plugins to define new commands. Because it is written in pure python, it is also very cross-platform, portable, easy to extend, and has no closed-source dependencies.

I've been meaning to write a "compiler" for amino acids, RNA and DNA for a while, and I have created a number of tools previously which I intend to embed in fastac as native functions, including PySplicer (which performs far more evidence-based and hopefully effective reverse translation than the stupid function I have included in fastac already), and DNAmespace, which provides a namespace-like interface to Genbank files at the genome-scale. With these tools and a bash-scripting like interface that can compile and use "libraries" of parts, a proper DNA "coding environment" will hopefully emerge.

My intention is that these tools, or their (better designed) successors, will help fill the gap between copy/paste and GUI-patina which has remained fairly untouched so far.

## What next?
Well, PySplicer needs a rewrite or a refactor for starters, and then I plan to integrate PySplicer properly into this "language" to provide a real reverse translation function that might actually work.

Then, I plan to integrate DNAmespace so that Genbank files can be used as Libraries, not just fasta files. This will allow richer access to genetic data such as by permitting "imports" of genes from entire bacterial genomes trivially.

I also plan to revamp the "macros" system a little so that macros are defined in a separate, more easily extended, file, or perhaps are imported from all files in a directory, allowing drag/drop addition of plugins.

I am happy to accept offers to extend this system provided they pass all the test cases and add useful functions. I am not interested in "coding style" or "PEP8" submissions, as I'm quite comfortable remembering the code as I have written it, thank you. :)
