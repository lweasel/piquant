import os
import os.path
import stat
import textwrap


class _Writer:
    def __init__(self):
        self.lines = []

    def _add_line(self, line_string):
        self.lines.append(line_string)

    def write_to_file(self, directory, filename):
        with open(os.path.join(directory, filename), "w") as f:
            f.write("\n".join(self.lines) + '\n')


class FluxSimulatorParamsWriter(_Writer):
    def __init__(self, vars_dict):
        _Writer.__init__(self)

        for name, value in vars_dict.items():
            self._add_line("{n} {v}".format(n=name, v=value))


class _BashSection:
    def __init__(self, writer):
        self.writer = writer

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.writer._add_line("")


class _BashBlock:
    def __init__(self, writer, start_prefix, start_suffix,
                 details, end, predeindent):
        self.writer = writer
        self.start_prefix = start_prefix
        self.start_suffix = start_suffix
        self.details = details
        self.end = end
        self.predeindent = predeindent

    def __enter__(self):
        self.writer.add_line("{p}{d}{s}".format(
            p=self.start_prefix, d=self.details, s=self.start_suffix))
        self.writer.indent()
        return self

    def __exit__(self, type, value, traceback):
        if self.predeindent:
            self.writer.deindent()
        self.writer.add_line(self.end)
        if not self.predeindent:
            self.writer.deindent()


class BashScriptWriter(_Writer):
    INDENT = '    '

    def __init__(self):
        _Writer.__init__(self)

        self.indent_level = 0
        self.block_ends = []

        with self.section():
            self.add_line("#!/bin/bash")
        with self.section():
            self.add_line("set -o nounset")
            self.add_line("set -o errexit")

    def indent(self):
        self.indent_level += 1

    def deindent(self):
        self.indent_level -= 1

    def add_line(self, line_string):
        line_string = BashScriptWriter.INDENT * self.indent_level + line_string
        _Writer._add_line(self, line_string)

    def section(self):
        return _BashSection(self)

    def block(self, block_start_prefix, block_start_suffix,
              block_end, details, predeindent=True):

        return _BashBlock(self, block_start_prefix, block_start_suffix,
                          details, block_end, predeindent)

    def if_block(self, details):
        return self.block("if [ ", " ]; then", "fi", details)

    def while_block(self, details):
        return self.block("while ", "; do", "done", details)

    def case_block(self, details):
        return self.block("case ", " in", "esac", details)

    def case_option_block(self, option):
        return self.block("", ")", ";;", option, predeindent=False)

    def add_comment(self, comment):
        lines = textwrap.wrap(
            comment, initial_indent="# ", subsequent_indent="# ",
            width=75 - self.indent_level * len(BashScriptWriter.INDENT))
        for line in lines:
            self.add_line(line)

    def add_echo(self, text=""):
        self.add_line("echo {t}".format(t=text))

    def add_pipe(self, pipe_commands):
        self.add_line(" | ".join(pipe_commands))

    def set_variable(self, variable, value):
        self.add_line("{var}={val}".format(var=variable, val=value))

    def write_to_file(self, directory, filename):
        _Writer.write_to_file(self, directory, filename)

        path = os.path.abspath(os.path.join(directory, filename))
        os.chmod(path,
                 stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
                 stat.S_IRGRP | stat.S_IROTH)
