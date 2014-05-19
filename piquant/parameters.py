import itertools
import ordutils.options as opt
import piquant_options as popt

QUANT_METHOD = "quant_method"
READ_DEPTH = "read_depth"
READ_LENGTH = "read_length"
PAIRED_END = "paired_end"
ERRORS = "errors"
BIAS = "bias"
POLYA_TAIL = "polya"


class _Parameter:
    def __init__(self, name, title,
                 option_name, option_validator,
                 is_numeric=False, value_namer=None, file_namer=None):
        self.name = name
        self.title = title
        self.option_name = option_name
        self.option_validator = option_validator
        self.is_numeric = is_numeric
        self.value_namer = value_namer if value_namer else lambda x: x
        self.file_namer = file_namer if file_namer else self.value_namer

    def get_value_name(self, value):
        return self.value_namer(value)

    def get_file_name_part(self, value):
        return self.file_namer(value)

_PARAMETERS = []

# TODO: get rid of isinstance() ASAP
_PARAMETERS.append(_Parameter(
    QUANT_METHOD, "Method", "--quant-method",
    popt.check_quantification_method,
    value_namer=lambda x: x if isinstance(x, basestring) else x.get_name()))

_PARAMETERS.append(_Parameter(
    READ_DEPTH, "Read depth", "--read-depth", int,
    is_numeric=True,
    value_namer=lambda x: "{d}x".format(d=x)))

_PARAMETERS.append(_Parameter(
    READ_LENGTH, "Read length", "--read-length", int,
    is_numeric=True,
    value_namer=lambda x: "{l}b".format(l=x)))

_PARAMETERS.append(_Parameter(
    PAIRED_END, "End type", "--paired-end",
    opt.check_boolean_value,
    value_namer=lambda x: "paired-end" if x else "single-end",
    file_namer=lambda x: "pe" if x else "se"))

_PARAMETERS.append(_Parameter(
    ERRORS, "Error type", "--error",
    opt.check_boolean_value,
    value_namer=lambda x: "with errors" if x else "no errors",
    file_namer=lambda x: "errors" if x else "no_errors"))

_PARAMETERS.append(_Parameter(
    BIAS, "Bias", "--bias",
    opt.check_boolean_value,
    value_namer=lambda x: "with bias" if x else "no bias",
    file_namer=lambda x: "bias" if x else "no_bias"))

_PARAMETERS.append(_Parameter(
    POLYA_TAIL, "PolyA tail", "--polya",
    opt.check_boolean_value,
    value_namer=lambda x: "with polyA tail" if x else "no polyA tail",
    file_namer=lambda x: "polya" if x else "no_polya"))


def get_parameters():
    return set(_PARAMETERS)


def validate_command_line_parameter_sets(options):
    param_vals = {}
    for param in get_parameters():
        if param.option_name in options:
            param_vals[param.name] = set(opt.validate_list_option(
                options[param.option_name], param.option_validator))
    return param_vals


def get_file_name(**params):
    elements = []
    for param in _PARAMETERS:
        if param.name in params:
            value = params[param.name]
            elements.append(param.get_file_name_part(value))
    return "_".join(elements)


def execute_for_param_sets(callables, **params_values):
    param_names = params_values.keys()
    for to_call in callables:
        for param_set in itertools.product(*params_values.values()):
            param_map = dict(zip(param_names, param_set))
            to_call(**param_map)
