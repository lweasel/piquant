import itertools
import options as opt
import quantifiers
import schema
import six

_PARAMETERS = []
_RUN_PARAMETERS = []


class _Parameter():
    def __init__(self, name, title, option_name, option_validator,
                 is_numeric=False, value_namer=None, file_namer=None,
                 run_parameter=True):

        self.name = name
        self.title = title
        self.option_name = option_name
        self.option_validator = option_validator
        self.run_parameter = run_parameter
        self.is_numeric = is_numeric
        self.value_namer = value_namer if value_namer else lambda x: x
        self.file_namer = file_namer if file_namer else self.value_namer

        _PARAMETERS.append(self)
        if run_parameter:
            _RUN_PARAMETERS.append(self)

    def get_value_name(self, value):
        return self.value_namer(value)

    def get_file_name_part(self, value):
        return self.file_namer(value)


TRANSCRIPT_GTF = _Parameter(
    "transcript_gtf", "Transcript GTF file", "--transcript-gtf",
    lambda x: opt.validate_file_option(
        x, "Transcript GTF file does not exist"),
    run_parameter=False)

GENOME_FASTA_DIR = _Parameter(
    "genome_fasta", "Genome FASTA directory", "--genome-fasta",
    lambda x: opt.validate_dir_option(
        x, "Genome FASTA directory does not exist"),
    run_parameter=False)

NUM_FRAGMENTS = _Parameter(
    "num_fragments", "Number of fragments", "--num-fragments",
    lambda x: opt.validate_int_option(
        x, "Number of fragments must be a positive integer", nonneg=True),
    run_parameter=False)

# TODO: get rid of isinstance() ASAP
is_string = lambda x: isinstance(x, six.string_types)

QUANT_METHOD = _Parameter(
    "quant_method", "Quantifier", "--quant-method",
    lambda x: quantifiers.get_quantification_methods()[x],
    value_namer=lambda x: x if is_string(x) else x.get_name())

READ_DEPTH = _Parameter(
    "read_depth", "Read depth", "--read-depth", int,
    is_numeric=True,
    value_namer=lambda x: "{d}x".format(d=x))

READ_LENGTH = _Parameter(
    "read_length", "Read length", "--read-length", int,
    is_numeric=True,
    value_namer=lambda x: "{l}b".format(l=x))

PAIRED_END = _Parameter(
    "paired_end", "End type", "--paired-end",
    opt.check_boolean_value,
    value_namer=lambda x: "paired-end" if x else "single-end",
    file_namer=lambda x: "pe" if x else "se")

ERRORS = _Parameter(
    "errors", "Error type", "--error",
    opt.check_boolean_value,
    value_namer=lambda x: "with errors" if x else "no errors",
    file_namer=lambda x: "errors" if x else "no_errors")

BIAS = _Parameter(
    "bias", "Bias", "--bias",
    opt.check_boolean_value,
    value_namer=lambda x: "with bias" if x else "no bias",
    file_namer=lambda x: "bias" if x else "no_bias")


def get_run_parameters():
    return set(_RUN_PARAMETERS)


# TODO: add tests for params_file and ignore_params options
def validate_command_line_parameter_sets(
        params_file, cl_options, ignore_params=[]):

    file_param_vals = {}
    if params_file:
        with open(params_file) as f:
            file_param_vals = {param: vals for param, vals in
                               [line.strip().split() for line in f]}

    param_vals = {}
    for param in _PARAMETERS:
        if param in ignore_params:
            continue

        for values_dict in [file_param_vals, cl_options]:
            if param.option_name in values_dict and \
                    values_dict[param.option_name] is not None:
                validated_vals = opt.validate_options_list(
                    values_dict[param.option_name], param.option_validator)
                param_vals[param.name] = set(validated_vals) \
                    if param.run_parameter else validated_vals[0]

        if param.name not in param_vals:
            raise schema.SchemaError(
                None, param.title + " parameter values must be specified.")

    return param_vals


def get_file_name(**params):
    elements = []
    for param in _RUN_PARAMETERS:
        if param.name in params:
            value = params[param.name]
            elements.append(param.get_file_name_part(value))
    return "_".join(elements)


def execute_for_param_sets(callables, **params_values):
    all_run_param_names = [p.name for p in _RUN_PARAMETERS]

    run_param_values = {}
    non_run_param_values = {}
    for param, values in params_values.items():
        if param in all_run_param_names:
            run_param_values[param] = values
        else:
            non_run_param_values[param] = values

    run_param_names = run_param_values.keys()
    for to_call in callables:
        for param_set in itertools.product(*run_param_values.values()):
            param_map = dict(zip(run_param_names, param_set))
            param_map.update(non_run_param_values)
            to_call(**param_map)
