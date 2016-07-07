# Module ValueExtractor

import collections
import inspect
import itertools
import re

####################################################################
# Metaclass for ordering items
####################################################################
class OrderedClass(type):
    """metaclass for preserving ordering of declarations"""
    @classmethod
    def __prepare__(metacls, name, bases):
        return collections.OrderedDict()
    def __new__(self, name, bases, classdict):
        classdict['__ordered__'] = [key for key in classdict.keys()
                if key not in ('__module__', '__qualname__')]
        return type.__new__(self, name, bases, classdict)

####################################################################
# Field types and their bases
####################################################################

class FieldBase:
    """Base class type for a field.

    Implementors must define the following:
    *  _default_format - the string used to format a single value. Must be
        at most 10 characters wide.
    """
    _counter = itertools.count()
    def __init__(self, format_string=None, set_context=False):
        """On printing, format_string is used to format if not None (otherwise
        use _default_format)."""
        self.order = next(FieldBase._counter)
        self.format_string = format_string if format_string else self._default_format
        self.set_context = set_context
        self.max_num_items = 1

class Expression(FieldBase):
    """Field that applies a python function to calculate a derived value.
    The function func() is examined and any arguments with variable names
    are replaced by those values."""

    _default_format = '%10.2f'

    def __init__(self, func, format_string=None, set_context=False):
        FieldBase.__init__(self, format_string, set_context)
        self.func = func

    def evaluate(self, context, results):
        sig = inspect.signature(self.func)
        args = []
        for name, param in sig.parameters.items():
            if name in results:
                if len(results[name])==1:
                    args.append(results[name][0])
                else:
                    args.append(results[name])
            elif name in context:
                if len(context[name])==1:
                    args.append(context[name][0])
                else:
                    args.append(context[name])
            else:
                if param.default == inspect.Parameter.empty: return None
                args.append(param.default)
        return [self.func(*args)]

class ExtractorMixin(FieldBase):
    """Provides basic search for a pattern and extract value functionality.

    Implementors must define the following:
    *  _default_format - the string used to format a single value. Must be
        at most 10 characters wide.
    * _default_add_pattern - the pattern added to the end of matchpattern if
        no groups are captured.
    """
    def __init__(self, matchpattern, format_string=None, set_context=False):
        """Field will match against supplied matchpattern. On printing,
        format_string is used to format if not None (otherwise use
        _default_format)."""
        FieldBase.__init__(self, format_string, set_context)
        self.max_num_items = 0 # Maximum number of items parsed
        self.matchpattern = matchpattern
        self.matchre = re.compile(self.matchpattern)
        if self.matchre.groups > 1:
            raise ValueError('matchpattern must have at most one group')
        if self.matchre.groups == 0:
            self.matchpattern = self.matchpattern + self._default_add_pattern
            self.matchre = re.compile(self.matchpattern)

# Represents a single number to be extracted
class SingleNumberMixin(ExtractorMixin):
    """Field expects a single number"""
    def parse(self, line):
        """Parse line to detect this pattern"""
        match = self.matchre.search(line)
        if match: self.max_num_items = 1
        return match.group(1) if match else None

class MultipleNumberMixin(ExtractorMixin):
    """Field expects multiple numbers"""
    def __init__(self, matchpattern, splitpattern=r'\s+', format_string=None,
            set_context=False):
        """In addition to searching on matchpattern, use splitpattern to split
        the group it identifies into multiple values."""
        ExtractorMixin.__init__(self, matchpattern, format_string, set_context)
        self.splitpattern = splitpattern
        self.splitre = re.compile(self.splitpattern)
    
    def parse(self, line):
        """Parse line to detect this pattern"""
        match = self.matchre.search(line)
        if not match: return None
        v = self.splitre.split(match.group(1))
        self.max_num_items = max(self.max_num_items, len(v))
        return v

class SingleInt(SingleNumberMixin):
    """Parse a single integer"""
    _default_format = '%10d'
    _default_add_pattern = r'[^0-9]*([0-9]*)'
    def parse(self, line):
        v = SingleNumberMixin.parse(self,line)
        return [int(v)] if v else None

class SingleFloat(SingleNumberMixin):
    """Parse a single float"""
    _default_format = '%10.2f'
    _default_add_pattern = r'[^0-9]*([0-9.E+-]*)'
    def parse(self, line):
        v = SingleNumberMixin.parse(self,line)
        return [float(v)] if v else None

class MultipleInt(MultipleNumberMixin):
    """Parse multiple integers"""
    _default_format = '%10d'
    _default_add_pattern = r'[^0-9]*([0-9 ]*)'
    def parse(self, line):
        v = MultipleNumberMixin.parse(self,line)
        return [int(x) for x in v] if v else None

class MultipleFloat(MultipleNumberMixin):
    """Parse multiple floats"""
    _default_format = '%10.2f'
    _default_add_pattern = r'[^0-9]*([0-9.E+ -]*)'
    def parse(self, line):
        v = MultipleNumberMixin.parse(self,line)
        return [float(x) for x in v if x] if v else None

####################################################################
# File parsing
####################################################################

def _isFileParser(cls):
    """Return true if cls is a FileParser."""
    return cls and inspect.isclass(cls) and issubclass(cls,FileParser)

class FileParser(metaclass=OrderedClass):
    """Base class for anything that parses a file using field extractors."""
    _counter = itertools.count()
    def __init__(self):
        self.order = next(FileParser._counter)
        self.fields = {x: getattr(self,x) for x in self.__ordered__ if isinstance(getattr(self,x),FieldBase)}
        self.field_order = list(self.fields.keys())
        self.field_order.sort(key=lambda f: self.fields[f].order)

    def parse(self, f, context):
        result = {}

        for line in f:
            for fname in self.field_order:
                ex = self.fields[fname]
                if not isinstance(ex, ExtractorMixin): continue
                val = ex.parse(line)
                if val:
                    result[fname] = val
                    if ex.set_context: context[fname] = val
        for fname in self.field_order:
            expr = self.fields[fname]
            if not isinstance(expr, Expression): continue
            val = expr.evaluate(context, result)
            if val:
                result[fname] = val
                if expr.set_context: context[fname] = val
        return result

    def has_field(self, field):
        """Returns true if field is a valid field name"""
        return (field in self.fields)

    def get_max_num_items(self, field):
        """Returns maximum number of items in field that have been parsed"""
        if field in self.fields:
            return self.fields[field].max_num_items
        else:
            return 0

    def get_format_string(self, field):
        return self.fields[field].format_string

####################################################################
# Series parsing
####################################################################

class SeriesParser(metaclass=OrderedClass):
    """Describes an application of one or more FileParsers to a collection of
    problems."""
    def __init__(self, rundir, problems):
        """Parses list of problems with files located in rundir."""
        self.problems = problems
        self.results = {}
        self.parsers = {x: getattr(self,x)() for x in self.__ordered__ if _isFileParser(getattr(self,x))}
        self.ext_order = list(self.parsers.keys())
        self.ext_order.sort(key=lambda ext: self.parsers[ext].order)
        for prblm in self.problems:
            r = {}
            context = {}
            filebase = re.sub(r'/', r'__', prblm)
            for pname in self.ext_order:
                p = self.parsers[pname]
                with open(rundir+'/'+filebase+p.extension) as f:
                    r[pname] = p.parse(f, context)
            self.results[prblm] = r

    def print(self, fields=None):
        """Prints results in a text table."""
        if not fields:
            print("Auto field extraction to be implemented...")
            return

        # Print header
        print("# Problem           ", end="")
        for f in fields:
            prblm = self.problems[0]
            nitem = 0
            for lbl, p in self.parsers.items():
                nitem += p.get_max_num_items(f)
            print(('%'+str(nitem*10)+'s') % f, end="")
        print()
        print("                    ", end="")
        for f in fields:
            prblm = self.problems[0]
            for ext in self.ext_order:
                itr = 0
                nitem = self.parsers[ext].get_max_num_items(f)
                if nitem == 0: continue
                idx = (nitem > 1)
                for x in range(nitem):
                    if idx:
                        print('%10s' % (ext+'['+str(itr)+']'), end="")
                    else:
                        print('%10s' % ext, end="")
                    itr = itr + 1
        print()


        # Print values
        for prblm in self.problems:
            print("%-20s" % prblm, end="")
            for f in fields:
                for ext in self.ext_order:
                    if not self.parsers[ext].has_field(f): continue
                    results = self.results[prblm][ext]
                    if f in results:
                        format_string = \
                            self.parsers[ext].get_format_string(f)
                        for x in results[f]:
                            print(format_string % x, end="")
                    else:
                        print("%10s" % '-', end="")
            print()
        #print(self.results)
