import os
import os.path
import stat


class _Writer:
    def __init__(self, initial_lines=[]):
        self.lines = list(initial_lines)

    def add_line(self, line_string):
        self.lines.append(line_string)

    def write_to_file(self, directory, filename):
        filename = directory + os.path.sep + filename
        with open(filename, "w") as f:
            f.write("\n".join(self.lines) + '\n')


class FluxSimulatorParamsWriter(_Writer):
    def __init__(self, vars_dict):
        lines = ["{n} {v}".format(n=name, v=value)
                 for name, value in vars_dict.items()]
        _Writer.__init__(self, initial_lines=lines)


class BashScriptWriter(_Writer):
    def __init__(self, vars_dict):
        self.vars_dict = vars_dict
        self.indent_level = 0
        _Writer.__init__(self, initial_lines=["#!/bin/bash", ""])

    def indent(self):
        self.indent_level += 1

    def deindent(self):
        self.indent_level -= 1

    def add_line(self, line_string):
        line_string = '\t' * self.indent_level + \
            line_string.format(**self.vars_dict)
        _Writer.add_line(self, line_string)

    def add_break(self):
        _Writer.add_line(self, "")

    def write_to_file(self, filename):
        _Writer.write_to_file(self, filename)
        path = os.path.abspath(filename)
        os.chmod(path,
                 stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
                 stat.S_IRGRP | stat.S_IROTH)
