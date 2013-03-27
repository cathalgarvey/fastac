#!/usr/bin/env python3
'''A simple "compiler" for commented fasta.'''
import argparse
import shlex
import json
import collections
from fastac import sequtils

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

def get_lib_var(string):
    '''Returns "lib" and "varname" for a given string of either "varname" or
    "lib.varname" form; lib defaults to None.
    String will often take the form "foo.fasta.some fasta block", in which case
    this function will return ("foo.fasta", "some fasta block"), as splitting
    commences from the right side of the string and stops after the first period
    it encounters.'''
    lib, varname = None, None
    if "." in string: lib, varname = string.rsplit(".",1)
    else: varname = string
    return lib, varname

Macros = {
# Place macro functions in this dictionary. They should accept a list of arguments
# as given by shlex.split: one convenient way to handle this is to define an
# argparse.ArgumentParser instance for each function and have the function call
# this ArgumentParser's parse_args() method on each invocation's list of arguments.
# Additionally each function should accept a dictionary that represents its immediate
# environment; at time of writing, this dict will contain previously parsed lines
# in the current block, and the Parser instance doing the parsing, allowing direct
# manipulation of Parser/Namespace data.
}

# imported_libs contains Parsers used to parse referenced "libraries", but not the
# current parser. Parsers are keyed by the name of their library file.
imported_libs = {}
def get_library(libname):
    if libname not in imported_libs:
        lib = FastaCompiler(Macros)
        lib.compile_file(libname)
        imported_libs[libname] = lib
    return imported_libs[libname]

def include(args, env_dict):
    # Using argparse allows flexible use of the argument list with optional args
    # etc, and suits the use of shlex.split() as it mimics a bash-like interface.
    # There is little point including "fluff" like description, help, etc,
    # as these macros will not be able to return help to the user!
    ArgP = argparse.ArgumentParser()
    ArgP.add_argument("block_name")
    ArgP.add_argument("--lib")
    args = ArgP.parse_args(args)
    # Below: if --lib is passed, then blockname is used exactly as given.
    # Otherwise, blockname is checked for a period character, which implies an
    # import, and is split by the final period character to give lib and blockname.
    # This means only single-depth imports are possible, so no lib1.lib2.fooblock,
    # but this is necessary as library files are specified by filename, and will
    # usually be called whatever.fasta, so imports are "whatever.fasta.someblock".
    if args.lib != None:
        libname, blockname = args.lib, args.block_name
    else:
        libname, blockname = get_lib_var(args.block_name)
    if libname:
        # If not already imported, import a multifasta "library" and use that as "lib".
        lib = get_library(libname)
    else:
        # Use current FastaCompiler object, passed as "namespace".
        lib = env_dict['namespace']
    return lib.get_block_sequence(blockname)
Macros['include'] = include

def _peer_call_include(block_name, lib_name, env_dict):
    'An internal shorthand for calling include from other "macros" despite pre-parsed args.'
    call_args = ["--lib", lib_name, block_name] if lib_name else [block_name]
    return Macros['include'](call_args, env_dict)
Macros['_peer_call_include'] = _peer_call_include

def complement(args, env_dict):
    ArgP = argparse.ArgumentParser()
    ArgP.add_argument("block_name")
    ArgP.add_argument("--lib")
    args = ArgP.parse_args(args)
    # This demonstrates trans-macro calls, but also the awkwardness of doing so
    seq = Macros['_peer_call_include'](args.block_name, args.lib, env_dict)
    seq = sequtils.get_complement(seq)
    return seq.lower()
Macros['complement'] = complement

def translate(args, env_dict):
    ArgP = argparse.ArgumentParser()
    ArgP.add_argument("block_name")
    ArgP.add_argument("--lib")
    ArgP.add_argument("--table", default="table1")
    args = ArgP.parse_args(args)
    seq = Macros['_peer_call_include'](args.block_name, args.lib, env_dict)
    aminoseq = sequtils.translate(seq, args.table)
    return aminoseq
Macros['translate'] = translate

def dumb_backtranslate(args, env_dict):
    ArgP = argparse.ArgumentParser()
    ArgP.add_argument("block_name")
    ArgP.add_argument("--lib")
    ArgP.add_argument("--table", default="table1")
    args = ArgP.parse_args(args)
    seq = Macros['_peer_call_include'](args.block_name, args.lib, env_dict)
    rtr_seq = sequtils.dumb_backtranslate(seq, args.table)
    return rtr_seq
Macros['dumb_backtranslate'] = dumb_backtranslate

def mutate(args, env_dict):
    '''Returns specified block with a single-point substitution.
    Usage: $mutate [--lib libfile] block position substitution'''
    ArgP = argparse.ArgumentParser()
    ArgP.add_argument("block_name")
    ArgP.add_argument("--lib")
    ArgP.add_argument("position", type=int)
    ArgP.add_argument("substitution", type=str)
    args = ArgP.parse_args(args)
    seq = Macros['_peer_call_include'](args.block_name, args.lib, env_dict)
    nseq = seq[:args.position-1] + args.substitution + seq[args.position:]
    return nseq.lower()
Macros['mutate'] = mutate

def def_template(args, env_dict):
    '''Registers the foregoing parsed_lines as a new template in the templates dictionary.
    Templates take the form of positional or named python format strings: positional
    arguments are parsed into a list and unpacked, while named arguments are searched for
    in the local namespace (the namespace is simply unpacked into the format arguments).
    Support for cross-library templating would be nice.
    Format:
    # No title block, or the template will also be registered in the "compiled" namespace.
    ; Stuff to precede first "argument":
    gcattgactagatc{0}
    ; Named import from current compiled namespace:
    {priorblock1}
    ; Another argument with more stuff to add after
    {1}cccaatctggtgctgtgt
    $def_template foo_template'''
    ArgP = argparse.ArgumentParser()
    ArgP.add_argument("templatename")
    args = ArgP.parse_args(args)
    env_dict['namespace'].templates[args.templatename] = ''.join(env_dict['current_lines'])
Macros['def_template'] = def_template

def use_template(args, env_dict):
    '''Should accept a variable number of string arguments which refer to blocks
    by name, then fetch those blocks' sequences as strings. These are used as
    positional arguments to the string format method. The namespace dict is also
    unpacked into the format method call, so format strings can embed local blocks
    by name as part of the template.'''
    # NB: Added a __str__ method to FastaBlock objects to allow them to be used
    #  directly in string format method, so can now unpack
    #  env_dict['namespace'].namespace directly with "**".
    ArgP = argparse.ArgumentParser()
    ArgP.add_argument("templatename")
    ArgP.add_argument("argblocks", nargs="+") # Result is a list of all free args.
    ArgP.add_argument("-r", "--raw", action="store_true")
    args = ArgP.parse_args(args)
    # Should write a getter for this so it can parse "foo.bar" or "foo:bar"
    # to get templates from libs.
    templatelib, templatename = get_lib_var(args.templatename)
    if templatelib:
        lib = get_library(templatelib)
    else:
        lib = env_dict['namespace']
    template = lib.templates[templatename]
    positional_seqs = []
    for blockname in args.argblocks:
        if args.raw:
            # In this mode, "blockname" is actually a literal string to use!
            positional_seqs.append(blockname)
        else:
            # The default mode; fetch blocks by name.
            # Now that include natively supports implicit imports like libname.fasta.blockname,
            # the use of library blocks to *fill* a template is supported, also.
            positional_seqs.append(Macros['_peer_call_include'](blockname, None, env_dict))
    return template.format(*positional_seqs, **lib.namespace)
Macros['use_template'] = use_template

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
        self.sequence = sequence.lower()
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

    def __str__(self):
        'Returns sequence only; this allows objects to be unpacked into string format method.'
        return self.sequence

class FastaCompiler(object):
    '''Contains methods for compiling blocks, multifasta files.
    Also acts as a scope for precompiled blocks.'''
    def __init__(self, macros={}, linewrap=50, lettercase="lower", namespace = {}, templates = {}):
        self.macros = macros
        self.linewrap = linewrap
        self.lettercase = lettercase
        self.namespace = collections.OrderedDict(namespace)
        self.templates = templates

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
                #raise E
                errmsg = "Error compiling block with first line "+Block.splitlines()[0]+":\n\t"+str(E)
                raise FastaCompileError(errmsg)

    def get_block(self, title):
        if title not in self.namespace:
            errmsg="Could not find {} - Current namespace: {}".format(title, str(self.namespace))
            raise ValueError(errmsg)
        return self.namespace[title]

    def get_block_sequence(self, title):
        return self.get_block(title).sequence

    def do_macro(self, macroline, current_lines):
        '''Is passed the macro call line and all lines already parsed.
        As macros are passed this and the Parser object itself, macros can
        independently define actions to take directly on the namespace or Parser.'''
        macroline = shlex.split(macroline.strip()[1:])
        # Passing a dict of environment stuff allows extension of environment
        # variables or objects passed to macros/functions without having to
        # rewrite them all again..
        environment = {"current_lines":current_lines,
                       "namespace":self}
        result = ''
        if macroline[0] in self.macros:
            # Macros should be passed the Compiler or Namespace object:
            result = self.macros[macroline[0]](macroline[1:], environment)
        else:
            errmsg = "Could not find macro/function named '{0}'".format(macroline[0])
            raise FastaCompileError(errmsg)
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
                result = self.do_macro(line, lines)
                # Not all macros may return results, but if they do, it's to be
                # included in current block.
                if result: lines.append(result)
            else:
                lines.append(line)
        # Make & Return / Register FastaObj if a title was found.
        # in the absence of a title, just drop everything and carry on, was
        # probably an anonymous block for macro calls or template definitions.
        if title:
            FastaObj = FastaBlock(title, ''.join(lines), meta)
            # Finally:
            if returnblock: return FastaObj
            else: self.namespace[FastaObj.title] = FastaObj

    def as_multifasta(self, preserve_meta=True, print_all=False):
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
            if "private" in self.namespace[subblock].meta and\
             self.namespace[subblock].meta["private"] and\
             not print_all: continue
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
    output = LocalCompiler.as_multifasta(Args.plain, Args.print_all)
    if Args.last:
        # Split to get only the last block.
        output = output.split("\n\n").pop()
    if Args.output:
        with open(Args.output, 'w') as OutFile:
            OutFile.write(output)
    else:
        print(output)

if __name__ == "__main__":
    ArgP = argparse.ArgumentParser(description="A simple 'compiler' for commented fasta.")
    ArgP.add_argument("fastafile", help="File to compile.")
    ArgP.add_argument("-o", "--output", help="Filename to save output to. Defaults to standard output.")
    ArgP.add_argument("-l", "--linelength", type=int, default=50,
                  help="Length to wrap sequence blocks around. Default is 50.")
    ArgP.add_argument("-c", "--case", default="lower",
                  help="Casing to present sequence in. Can be either 'lower' or 'upper'. Defaults to lower.")
    ArgP.add_argument("-p", "--plain", default=True, action="store_false",
                  help="Output plain FASTA without metadata in title line.")
    ArgP.add_argument("-P", "--print-all", default=False, action="store_true",
                  help="Prints all blocks, including those with the 'private' metatag.")
    ArgP.add_argument("-L", "--last", default=False, action="store_true",
                  help="Only output the last fasta block compiled in the main file.")
    main(ArgP.parse_args())
