import itertools
import ordutils.options as opt
import quantifiers
import schema


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
QUANT_METHOD = _Parameter(
    "quant_method", "Method", "--quant-method",
    lambda x: quantifiers.get_quantification_methods()[x],
    value_namer=lambda x: x if isinstance(x, basestring) else x.get_name())
_PARAMETERS.append(QUANT_METHOD)

READ_DEPTH = _Parameter(
    "read_depth", "Read depth", "--read-depth", int,
    is_numeric=True,
    value_namer=lambda x: "{d}x".format(d=x))
_PARAMETERS.append(READ_DEPTH)

READ_LENGTH = _Parameter(
    "read_length", "Read length", "--read-length", int,
    is_numeric=True,
    value_namer=lambda x: "{l}b".format(l=x))
_PARAMETERS.append(READ_LENGTH)

PAIRED_END = _Parameter(
    "paired_end", "End type", "--paired-end",
    opt.check_boolean_value,
    value_namer=lambda x: "paired-end" if x else "single-end",
    file_namer=lambda x: "pe" if x else "se")
_PARAMETERS.append(PAIRED_END)

ERRORS = _Parameter(
    "errors", "Error type", "--error",
    opt.check_boolean_value,
    value_namer=lambda x: "with errors" if x else "no errors",
    file_namer=lambda x: "errors" if x else "no_errors")
_PARAMETERS.append(ERRORS)

BIAS = _Parameter(
    "bias", "Bias", "--bias",
    opt.check_boolean_value,
    value_namer=lambda x: "with bias" if x else "no bias",
    file_namer=lambda x: "bias" if x else "no_bias")
_PARAMETERS.append(BIAS)


def get_parameters():
    return set(_PARAMETERS)


# TODO: add tests for params_file and ignore_params options
def validate_command_line_parameter_sets(
        params_file, cl_options, ignore_params=[]):

    file_param_vals = {}
    if params_file:
        with open(params_file) as f:
            file_param_vals = {param: vals for param, vals in
                               [line.strip().split() for line in f]}

    param_vals = {}
    for param in get_parameters():
        if param.name in ignore_params:
            continue

        for values_dict in [file_param_vals, cl_options]:
            if param.option_name in values_dict and \
                    values_dict[param.option_name] is not None:
                param_vals[param.name] = set(opt.validate_list_option(
                    values_dict[param.option_name], param.option_validator))

        if param.name not in param_vals:
            raise schema.SchemaError(
                None, param.title + " parameter values must be specified.")

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
