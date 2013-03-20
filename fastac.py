#!/usr/bin/env python3
'''A simple "compiler" for commented fasta.'''
import argparse
import shlex
import sequtils
import json
import collections

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
            # Spawn a new FastaCompiler to import target library file.
            lib = FastaCompiler(Macros)
            lib.compile_file(args.lib)
            # Remember this lib in case it's referred again.
            imported_libs[args.lib] = lib
        # If already imported, use preexisting FastaCompiler object from imported_libs.
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

class FastaBlock(object):
    FastaFormat = "> {0}\n{1}"
    def __init__(self, title, sequence, meta):
        self.title = title
        self.sequence = sequence
        self.meta = meta
        if "type" in self.meta:
            self.type = self.meta['type']
        else:
            self.type = sequtils.deduce_alphabet(self.sequence)
            self.meta['type'] = self.type

    @staticmethod
    def _chunks(l, n):
        for i in range(0, len(l), n): yield l[i:i+n]

    def as_dict(self):
        'Used for exporting whole namespaces as JSON by converting subunits to dicts.'
        return {"title":self.title, "sequence":self.sequence,
                           "meta":self.meta, "type":self.type}

    def as_json(self, indent=2):
        'Returns a json string defining the content of this object for easy export/import.'
        return json.dumps(self.as_dict(), indent=indent)

    def as_fasta(self, linewrap=50):
        'Exports as vanilla FASTA.'
        OutSeq = '\n'.join([x for x in self._chunks(self.sequence, linewrap)])
        return self.FastaFormat.format(self.title, OutSeq)

    def as_metafasta(self, linewrap=50):
        'Exports as valid FASTA, but preserves metadata as a (huge, ugly) title line extension.'
        Metatitle = '{} {}'.format(self.title, json.dumps(self.meta))
        OutSeq = '\n'.join([x for x in self._chunks(self.sequence, linewrap)])
        return self.FastaFormat.format(Metatitle, OutSeq)

class FastaCompiler(object):
    '''Contains methods for compiling blocks, multifasta files.
    Also acts as a scope for precompiled blocks.'''
    def __init__(self, macros={}, linewrap=50, lettercase="lower", namespace = {}):
        self.macros = macros
        self.linewrap = linewrap
        self.lettercase = lettercase
        self.namespace = collections.OrderedDict(namespace)
    def compile_file(self, filen):
        with open(filen) as InputFile:
            self.compile_multifasta(InputFile.read())

    def compile_multifasta(self, file_contents):
        for Block in file_contents.strip().split("\n\n"):
            Block = Block.strip() # Handles more than one blank line b/w blocks.
            try:
                self.compile_block(Block)
            except Exception as E:
                # For debug, just raise. Can later sort out common exceptions
                # and raise more informative errors or catch/ignore.
                raise E

    def get_block(self, title):
        if title not in self.namespace:
            errmsg="Could not find {} - Current namespace: {}".format(title, str(self.namespace))
            raise ValueError(errmsg)
        return self.namespace[title]

    def get_block_sequence(self, title):
        return self.get_block(title).sequence

    def do_macro(self, macroline):
        '''Handles macro calls. Can be made more useful/complex later by
        allowing macros to define "actions" to undertake with returned data,
        so that macros can return entire sequence blocks to be registered in
        the namespace directly rather than included in the current compile
        block, for example.'''
        macroline = shlex.split(macroline.strip()[1:])
        if macroline[0] in self.macros:
            # Macros should be passed the Compiler or Namespace object:
            result = self.macros[macroline[0]](macroline[1:], self)
        if result: return result

    @staticmethod
    def import_inline_meta(meta, titlemeta, mergeconflicts = True):
        '''Given a meta dict and JSON found in a title block, check if titlemeta
        is "valid" fastac json, and import any keys found.
        This tries to avoid clobbering in general, but titlemeta subkeys defining
        a dict which already exists in meta will be applied as a dict update, in
        which case overwriting may occur. Also, if both define a list with the
        same name (for example "comments"), then the titlemeta list will extend
        meta, in which case duplicates may occur.
        Finally, returns the sequence type, if defined, or an empty string.'''
        if not isinstance(titlemeta, dict): return ''
        for key in titlemeta:
            if key not in meta:
                meta[key] = titlemeta[key]
            else:
                # Merge conflicts (overwriting keys in case of two dicts)
                if isinstance(meta[key],list) and isinstance(titlemeta[key],list):
                    meta[key].extend(titlemeta[key])
                elif isinstance(meta[key],dict) and isinstance(titlemeta[key],dict):
                    meta[key].update(titlemeta[key])

    @staticmethod
    def handle_markup(line, pos, meta):
        '''Processes line, seeking a json comment or just using the whole line
        and pos to add a comment to "meta" dict.'''
        inline_meta, line = _getjson(line)
        if len(inline_meta) > 1:
            errmsg = ("Only one inline JSON item can be defined per "
                        "metadata line:\n\t"+line)
            raise FastaCompileError(errmsg)
        if "comment" in inline_meta:
            # Inline json comments take precedence and other line
            # contents are ignored entirely.
            # They take the form: {"comment":[4,55, "Promoter and RBS"]}
            if not isinstance(inline_meta['comment'], list) or\
                not isinstance(inline_meta['comment'][0], int) or\
                not isinstance(inline_meta['comment'][1], int) or\
                not isinstance(inline_meta['comment'][2], str):
                errmsg = ("JSON inline comments must be of form "
                        "[int, int, string]:\n\t"+line)
                raise FastaCompileError(errmsg)
            comment = inline_meta['comment']
        else:
            comment = [pos, pos, line.lstrip(";").lstrip()]
        meta['comments'].append(comment)

    def compile_block(self, block, returnblock=False):
        '''Compiles a FASTA sequence block, possibly with macros.
        If returnblock, then the resulting compiled FASTA object is returned.
        Otherwise, it is added to this FastaCompiler's namespace attribute.'''
        if not isinstance(block, str):raise ValueError("block must be a string")
        # As compile_fasta_block
        title = ''
        lines = []
        # Comment format is [int(start), int(finish), str(comment)], like [1,14,"Promoter"]
        meta = {"comments":[]}
        for line in block.splitlines():
            line = line.strip()
            if line[0] == ">":
                if title:
                    errmsg =("Title already defined for this block, but another"
                    " title line has been encountered:\n\t"+line)
                    raise FastaCompileError(errmsg)
                line = line.lstrip(">").lstrip()
                json_meta, title = _getjson(line)
                for json_object in json_meta:
                    # More than one may occur, although that would be dumb.
                    # import_title_meta extends or overwrites meta so no return
                    # is needed.
                    self.import_inline_meta(meta, json_object)
            elif line[0] == ";":
                # Positional sequence markup metadata.
                pos = len(''.join(lines)) + 1
                self.handle_markup(line, pos, meta)
            elif line[0] == "#":
                # Comments, don't keep.
                pass
            elif line[0] == "$":
                result = self.do_macro(line)
                # Not all macros may return results, but if they do, it's to be
                # included in current block.
                if result: lines.append(result)
            else:
                lines.append(line)
        # Make FastaObj:
        #print("Compiled:\n> {}\n{}".format(title, ''.join(lines)))
        FastaObj = FastaBlock(title, ''.join(lines), meta)
        # Finally:
        if returnblock: return FastaObj
        else: self.namespace[FastaObj.title] = FastaObj

    def as_multifasta(self, preserve_meta=True):
        '''Return namespace in order of compilation as a multi-fasta file.
        If preserve_meta is true, export as "metafasta", where metadata is
        kept in a json block in the title. This is ugly, but lossless and cross-
        compatible with other bioinfo tools, which will ignore the big title.'''
        Export_Blocks = []
        for subblock in self.namespace:
            # As self.namespace is an OrderedDict, this will export in the same
            # order as compilation occurred.
            # Could also achieve this effect by keeping a list of compiled titles
            # and iterating over the list to get things from the namespace dict
            # in order: might be worth checking if this is more efficient?
            if preserve_meta: Export_Blocks.append(self.namespace[subblock].as_metafasta())
            else: Export_Blocks.append(self.namespace[subblock].as_fasta())
        return '\n\n'.join(Export_Blocks)

    def as_json(self, indent=2):
        'This makes a more efficient library format so may be preferred in future.'
        jsonablenamespace = {}
        for FastaObject in self.namespace:
            jsonablenamespace[FastaObject.title] = FastaObject.as_dict()
        return json.dumps(jsonablenamespace, indent=indent)

def main(Args):
    'Expects an argparse parse_args namespace.'
    LocalCompiler = FastaCompiler(Macros, Args.linelength, Args.case)
    LocalCompiler.compile_file(Args.fastafile)
    if Args.output:
        with open(Args.output, 'w') as OutFile:
            OutFile.write(LocalCompiler.as_multifasta())
    else:
        print(LocalCompiler.as_multifasta())

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
