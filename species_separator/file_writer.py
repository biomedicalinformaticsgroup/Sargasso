import contextlib
import os
import os.path
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


class MakefileWriter(_Writer):
    INDENT = '\t'

    def __init__(self):
        _Writer.__init__(self)

        self.indent_level = 0

    def indent(self):
        self.indent_level += 1

    def deindent(self):
        self.indent_level -= 1

    def add_line(self, line_string):
        line_string = MakefileWriter.INDENT * self.indent_level + line_string
        _Writer._add_line(self, line_string)

    def add_blank_line(self):
        self.add_line("")

    @contextlib.contextmanager
    def target_definition(self, target, dependencies, raw_target=False, raw_dependencies=False):
        self.add_line("{tar}: {deps}".format(
            tar=self.variable_val(target, raw_target),
            deps=" ".join([self.variable_val(d, raw_dependencies) for d in dependencies])))
        self.indent()

        try:
            yield
        finally:
            self.deindent()
            self.add_blank_line()

    def add_comment(self, comment):
        lines = textwrap.wrap(
            comment, initial_indent="# ", subsequent_indent="# ",
            width=75 - self.indent_level * len(MakefileWriter.INDENT))
        for line in lines:
            self.add_line(line)

    def set_variable(self, variable, value):
        self.add_line("{var}={val}".format(var=variable, val=value))

    def variable_val(self, variable, raw=False):
        return variable if raw else "$({var})".format(var=variable)

    def make_target_directory(self, target):
        self.add_command("mkdir", ["-p", self.variable_val(target)])

    def remove_target_directory(self, target, raw_target=False):
        self.add_command("rm", ["-rf", self.variable_val(target, raw_target)])

    def add_command(self, command_name, options):
        line_elements = [command_name] + options
        self.add_line(" ".join(line_elements))
