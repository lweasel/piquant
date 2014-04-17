import os
import os.path
import stat
import textwrap


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
    INDENT = '    '

    def __init__(self, vars_dict):
        self.vars_dict = vars_dict
        self.indent_level = 0
        self.block_ends = []
        _Writer.__init__(self, initial_lines=["#!/bin/bash", ""])

    def indent(self):
        self.indent_level += 1

    def deindent(self):
        self.indent_level -= 1

    def add_line(self, line_string):
        line_string = BashScriptWriter.INDENT * self.indent_level + \
            line_string.format(**self.vars_dict)
        _Writer.add_line(self, line_string)

    def start_block(self, block_start_prefix, block_start_suffix,
                    block_end, details, predeindent=True):
        self.add_line("{p}{d}{s}".format(
            p=block_start_prefix, d=details, s=block_start_suffix))
        self.block_ends.append(block_end)
        self.indent()
        self.predeindent = predeindent

    def end_block(self):
        if self.predeindent:
            self.deindent()
        self.add_line(self.block_ends[-1])
        self.block_ends = self.block_ends[:-1]
        if not self.predeindent:
            self.deindent()
        self.predeindent = True

    def start_if(self, details):
        self.start_block("if [ ", " ]; then", "fi", details)

    def start_while(self, details):
        self.start_block("while ", "; do", "done", details)

    def start_case(self, details):
        self.start_block("case ", "in", "esac", details)

    def start_case_option(self, option):
        self.start_block("", ")", ";;", option, predeindent=False)

    def add_comment(self, comment):
        lines = textwrap.wrap(
            comment, initial_indent="# ", subsequent_indent="# ",
            width=75 - self.indent_level * len(BashScriptWriter.INDENT))
        for line in lines:
            self.add_line(line)

    def add_echo(self, text=""):
        self.add_line("echo {t}".format(t=text))

    def add_break(self):
        _Writer.add_line(self, "")

    def write_to_file(self, directory, filename):
        _Writer.write_to_file(self, directory, filename)
        path = os.path.abspath(directory + os.path.sep + filename)
        os.chmod(path,
                 stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
                 stat.S_IRGRP | stat.S_IROTH)
