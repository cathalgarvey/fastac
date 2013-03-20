#!/usr/bin/env python3
'''A simple "compiler" for commented fasta.'''
import argparse
import shlex
import sequtils
import json

# Handy functions:
def _chunks(l, n):
    "Yield successive n-sized chunks from l."
    for i in range(0, len(l), n): yield l[i:i+n]

def _case(string, lettercase="lower"):
    for x in [string, lettercase]:
        if not isinstance(x, str): raise TypeError("Args must be strings.")
    lettercase = lettercase.lower()
    if lettercase == "lower": return string.lower()
    elif lettercase == "upper": return string.upper()
    elif lettercase == "preserve": return string
    else: raise ValueError("Argument 'lettercase' can be either 'upper' or 'lower' or 'preserve'.")

def _getjson(string):
    '''Finds and returns a list of json objects found within a string.
    Silently ignores failed decodes, so substrings like {this} will not be
    decoded and returned, and will not trigger an Exception.'''
    blocks, strblocks, buf = [], [], []
    buffering = False
    for char in string:
        # Buffering approach, with buffer reset only taking place if the current
        # buffer can be decoded with the json.loads() function successfully.
        if char in "{[": buffering = True
        elif char in "]}":
            buf.append(char)
            try:
                newobj = json.loads(''.join(buf))
                strblocks.append(''.join(buf))
                blocks.append(newobj)
                buf, buffering = [], False
            except: continue
        if buffering: buf.append(char)
    # Now create and return a copy of the title without the json objects.
    for i in strblocks: string = string.replace(i, '').strip()
    return blocks, string

class MultiFasta(list):
    def __init__(self, *args, **nargs):
        list.__init__(self, *args, **nargs)
        self.blocks_by_title = {}
    def register(self, fasta_obj):
        # New item will be current length of list, as list indices start at 0
        newblock_index = len(self)
        self.append(fasta_obj)
        self.blocks_by_title[fasta_obj['Title']] = newblock_index
    def get_block(self, title):
        if title not in self.blocks_by_title:
            errmsg = "No precompiled fasta block could be found with this title: "+title
            errmsg += "\nPrecompiled blocks: "+str(self.blocks_by_title.keys())
            raise ValueError(errmsg)
        return self[self.blocks_by_title[title]]
    def get_block_sequence(self, title):
        block = self.get_block(title)
        return block['Sequence']
    def fasta_list(self, linewrap=50):
        FFormat = "> {0}\n{1}"
        Outblocks = []
        for i in self:
            OutSequence = '\n'.join(_chunks(i['Sequence'], linewrap))
            Outblocks.append(FFormat.format(i['Title'], OutSequence))
        return Outblocks

CompiledBlocks = MultiFasta()

Macros = {
# Place macro functions in this dictionary. They should accept a list of arguments
# as given by shlex.split: one convenient way to handle this is to define an
# argparse.ArgumentParser instance for each function and have the function call
# this ArgumentParser's parse_args() method on each invocation's list of arguments.
}

# Sample macro, either references a pre-compiled fasta object or imports a library
# as a new MultiFasta file and gets the fasta object from that.
imported_libs = {}
def include(args, currentnamespace):
    # Using argparse allows flexible use of the argument list with optional args
    # etc, and suits the use of shlex.split() perfectly as it mimics a bash-like
    # interface.
    # However, there is little point including "fluff" like description, help, etc,
    # as these macros are called as functions and will not be able to return help
    # to the user!
    ArgP = argparse.ArgumentParser()
    ArgP.add_argument("block_name")
    ArgP.add_argument("--lib")
    args = ArgP.parse_args(args)
    if args.lib:
        # If not already imported, import a multifasta "library" and use that as "lib".
        if args.lib not in imported_libs:
            with open(args.lib) as IncludeFile:
                libcontents = IncludeFile.read()
            lib = compile_multifasta(libcontents, MultiFasta(), objexport=True, macros=Macros)
            # Remember this lib in case it's referred again.
            imported_libs[args.lib] = lib
        # If already imported, use preexisting MultiFasta object from imported_libs.
        else: lib = imported_libs[args.lib]
    else:
        lib = currentnamespace
    return lib.get_block_sequence(args.block_name)
Macros['include'] = include

def complement(args, currentnamespace):
    ArgP = argparse.ArgumentParser()
    ArgP.add_argument("block_name")
    ArgP.add_argument("--lib")
    args = ArgP.parse_args(args)
    # This demonstrates trans-macro calls, but also the awkwardness of doing so with
    # optional arguments and argparse..
    inc_call = [args.block_name, "--lib", args.lib] if args.lib else [args.block_name]
    seq = Macros['include'](inc_call)
    seq = sequtils.get_complement(seq)
    return seq.lower()
Macros['complement'] = complement

def translate(args, currentnamespace):
    ArgP = argparse.ArgumentParser()
    ArgP.add_argument("block_name")
    ArgP.add_argument("--lib")
    ArgP.add_argument("--table", default="table1")
    args = ArgP.parse_args(args)
    inc_call = [args.block_name, "--lib", args.lib] if args.lib else [args.block_name]
    seq = Macros['include'](inc_call, currentnamespace)
    aminoseq = sequtils.translate(seq, args.table)
    return aminoseq
Macros['translate'] = translate

def dumb_backtranslate(args, currentnamespace):
    ArgP = argparse.ArgumentParser()
    ArgP.add_argument("block_name")
    ArgP.add_argument("--lib")
    ArgP.add_argument("--table", default="table1")
    args = ArgP.parse_args(args)
    inc_call = [args.block_name, "--lib", args.lib] if args.lib else [args.block_name]
    seq = Macros['include'](inc_call, currentnamespace)
    rtr_seq = sequtils.dumb_backtranslate(seq, args.table)
    return rtr_seq
Macros['dumb_backtranslate'] = dumb_backtranslate

class FastaError(Exception):
    'For errors deriving from bad fasta input.'
    pass

class FastaCompileError(FastaError):
    'For errors that arise as consequence of attempting to compile fasta.'
    pass

def compile_fasta_block(fasta_block, current_namespace, linewrap=50,
                        lettercase="lower", objexport=False, macros={}):
    FastaFormat = '> {Title}\n{Sequence}'
    fasta_block = fasta_block.strip()
    title = ''
    seqtype = ''
    lines = []
    meta = {'comments':[]}
    for line in fasta_block.splitlines():
        line = line.strip()
        if line[0] == ";":
            # Markup comment? Keep comment line + position as meta?
            # Get human-readable index of next character, which comment is assumed to
            # refer to.
            commented_index = len(''.join(lines)) + 1
            # list, not tuple, to avoid conflicts with JSON meta.
            comment = [commented_index, commented_index, line]
            meta['comments'].append(comment)
        elif line[0] == "#":
            # Unrecorded comment extension. Ignore line.
            continue
        elif line[0] == "$":
            # Command/Macro extension?
            macroline = shlex.split(line[1:])
            if macroline[0] in macros:
                result = macros[macroline[0]](macroline[1:], current_namespace)
                if result: lines.append(result)
        elif line[0] == ">":
            if title: raise FastaError("Found second title line in Fasta Block.")
            # Extracts any JSON metadata extensions in the title and returns the
            # title without them.
            jsoncontent, extracted_line = _getjson(line)
            if jsoncontent:
                # JSON content can contain metadata including sequence type.
                # JSON metadata that's compliant with this extension will contain
                # a "fastac" key with a version, starting at "1".
                # Compliant objects will then define a "type" and "meta" key,
                # respectively a string ("dna", "rna" or "amino") and a list of
                # meta keys to include.
                for json_item in jsoncontent:
                    if not isinstance(json_item, dict): continue
                    if "fastac" in json_item:
                        seqtype = json_item['type']
                        for item in json_item['meta']:
                            itemcontent = json_item['meta'][item]
                            #meta.append(item)
                            if item == "comments":
                                meta['comments'].extend(itemcontent)
                            else:
                                if item in meta.keys():
                                    errmsg= "Meta key '{}' already in meta dict but is being overwritten by title-meta.".format(item)
                                    raise ValueError(errmsg)
                                meta[item] = itemcontent
            title = extracted_line.lstrip(">").lstrip()
        else:
            lines.append(line)
    sequence = _case(''.join(lines), lettercase)
    if not title: raise FastaError("No title found for this block.")
    if not seqtype:
        # If not defined by meta, then guess sequence type by letter content.
        sequtils.deduce_alphabet(sequence)
    if not objexport:
        sequence = [x for x in _chunks(sequence, linewrap)]
        sequence = '\n'.join(sequence)
        return FastaFormat.format(Title=title, Sequence=sequence)
    else:
        return {'Title':title, 'Sequence':sequence, 'Meta':meta}

def compile_multifasta(file_contents, multifasta_container, linelength=50,
                       lettercase='lower', objexport=False, macros={}):
    'Splits by empty lines and passes each block to compile_fasta_block. Returns list of outputs.'
    # CompiledBlocks is defined in global scope, above, so that macros may use it.
    # Used in error reporting.
    current_block = 1
    for Block in file_contents.strip().split("\n\n"):
        try:
            # DoneBlock is passed current block, the current namespace (this
            # allows library-recursion), and args like linelength and lettercase.
            #
            DoneBlock = compile_fasta_block(Block, multifasta_container,
                        linelength, lettercase, objexport=True, macros=macros)
            multifasta_container.register(DoneBlock)
        except Exception as E:
            # For debug uncomment the following line:
            raise E
            # Count the newlines until the current block and call that the line number.
            error_line = file_contents[:file_contents.find(Block)].count("\n") #+ 1 # ?
            ErrorStr = "Error on line {}, FASTA block {}: ".format(error_line, current_block) + str(E)
            raise FastaCompileError(ErrorStr)
        current_block += 1
    if objexport:
        return multifasta_container
    else:
        return multifasta_container.fasta_list(linelength)

def main(Args):
    'Expects an argparse parse_args namespace.'
    with open(Args.fastafile) as InputFile:
        Compiled = compile_multifasta(InputFile.read(), CompiledBlocks,
                                            Args.linelength, Args.case, macros=Macros)
    if Args.output is not None:
        with open(Args.output, "w") as OutputFile:
            OutputFile.write('\n\n'.join(Compiled))
    else:
        print("\n\n".join(Compiled))

if __name__ == "__main__":
    ArgP = argparse.ArgumentParser(description="A simple 'compiler' for commented fasta.")
    ArgP.add_argument("fastafile", help="File to compile.")
    ArgP.add_argument("-o", "--output", help="Filename to save output to. Defaults to standard output.")
    ArgP.add_argument("-l", "--linelength",
                      help="Length to wrap sequence blocks around. Default is 50.",
                      type=int, default=50)
    ArgP.add_argument("-c", "--case", default="lower",
                  help="Casing to present sequence in. Can be either 'lower' or 'upper'. Defaults to lower.")
    main(ArgP.parse_args())
