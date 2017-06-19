import contextlib
import os
import os.path
import stat
import textwrap


@contextlib.contextmanager
def writing_to_file(writer_cls, directory, filename):
    writer = writer_cls()
    try:
        yield writer
    finally:
        writer.write_to_file(directory, filename)


class _Writer(object):
    def __init__(self):
        self.lines = []

    def _add_line(self, line_string):
        self.lines.append(line_string)

    def write_to_file(self, directory, filename):
        with open(os.path.join(directory, filename), "w") as output_file:
            output_file.write("\n".join(self.lines) + '\n')


class FluxSimulatorParamsWriter(_Writer):
    def __init__(self):
        _Writer.__init__(self)

    def add_vars(self, vars_dict):
        for name, value in vars_dict.items():
            self._add_line("{n} {v}".format(n=name, v=value))


class PolyesterSimulatorParamsWriter(_Writer):
    def __init__(self):
        _Writer.__init__(self)

    def add_vars(self, vars_dict):
        for name, value in vars_dict.items():
            self._add_line("{n} = {v}".format(n=name, v=value))


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

    @contextlib.contextmanager
    def section(self):
        try:
            yield
        finally:
            self._add_line("")

    def _adding_bash_block(
            self, start_prefix, start_suffix, end, details, predeindent=True):

        self.add_line("{p}{d}{s}".format(
            p=start_prefix, d=details, s=start_suffix))
        self.indent()

        try:
            yield
        finally:
            if predeindent:
                self.deindent()
            self.add_line(end)
            if not predeindent:
                self.deindent()

    @contextlib.contextmanager
    def if_block(self, test_command):
        return self._adding_bash_block("if [ ", " ]; then", "fi", test_command)

    @contextlib.contextmanager
    def while_block(self, details):
        return self._adding_bash_block("while ", "; do", "done", details)

    @contextlib.contextmanager
    def case_block(self, details):
        return self._adding_bash_block("case ", " in", "esac", details)

    @contextlib.contextmanager
    def case_option_block(self, option):
        return self._adding_bash_block(
            "", ")", ";;", option, predeindent=False)

    def add_comment(self, comment):
        lines = textwrap.wrap(
            comment, initial_indent="# ", subsequent_indent="# ",
            width=75 - self.indent_level * len(BashScriptWriter.INDENT))
        for line in lines:
            self.add_line(line)

    def add_echo(self, text=""):
        self.add_line("echo {t}".format(t=text))

    def add_pipe(self, *pipe_commands):
        self.add_line(" | ".join(pipe_commands))

    def set_variable(self, variable, value):
        self.add_line("{var}={val}".format(var=variable, val=value))

    def write_to_file(self, directory, filename):
        _Writer.write_to_file(self, directory, filename)

        path = os.path.abspath(os.path.join(directory, filename))
        os.chmod(path,
                 stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR |
                 stat.S_IRGRP | stat.S_IROTH)

class RScriptWriter(_Writer):
    def __init__(self):
        _Writer.__init__(self)

    def add_comment(self, comment):
        lines = textwrap.wrap(
                comment, initial_indent="# ",subsequent_indent="# ", width=70)
        for line in lines:
            self._add_line(line)

    def get_variable(self,value,isInitial):
        if isInitial:
            if isinstance(value, str):
                value = '"' + value + '"'
            else:
                if isinstance(value,bool):
                    value = "TRUE" if value else "FALSE"
        return value

    def add_variable(self,variable,value,isInitial):
        self._add_line("{var} = {val}".format(var=variable,val=self.get_variable(value,isInitial)))

    def call_function(self,function,par_dict):
        myline = function + "("
        for key,value in par_dict.items():
            if key.startswith("abbrev_"):
                myline += value + ",\n\t"
            else:
                myline += "{var} = {val},\n\t".format(var=key,val=self.get_variable(value,False))
        myline = myline[:-3] + ")"
        return myline

    def add_library(self,lib_list):
        for library in lib_list:
            self._add_line('library("{lib}")'.format(lib=library))

    def add_function_variable(self,variable,function,par_dict):
        value = self.call_function(function, par_dict)
        myline = "{var} = {val}".format(var=variable,val=self.get_variable(value,False))
        self._add_line(myline)

    def get_vector_value(self,val_list):
        #val_list = vals.split(',')
        myline = r'c('
        for value in val_list:
            myline += '"{val}",\n'.format(val=value)
        myline = myline[:-2] + r')'
        return myline

    def add_vector_variable(self,variable,val_list):
        self.add_variable(variable,self.get_vector_value(val_list),False)

    def get_dataframe_column(self,data_frame,column):
        return "{df}${col_name}".format(df=data_frame,col_name=column)

    def get_indices_specification(self,indices):
        myindices = r"["
        indices = indices.split(r",")
        for index in indices:
            myindices += index + r","
        myindices = myindices[:-1] + r"]"
        return myindices

    def get_expression(self,left,operator,right):
        return "{l} {o} {r}".format(l=left,o=operator,r=right)

    def get_long_expression(self,var_list,operator_list):
        my_expression = ""
        operator_list.append("")
        for i in range(len(var_list)):
            my_expression += operator_list[i-1] + var_list[i]
        return my_expression
