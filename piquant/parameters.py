class _Parameter:
    def __init__(self, name, title, is_numeric=False,
                 value_namer=None, file_namer=None):
        self.name = name
        self.title = title
        self.is_numeric = is_numeric
        self.value_namer = value_namer if value_namer else lambda x: x
        self.file_namer = file_namer if file_namer else self.value_namer

    def get_value_name(self, value):
        return self.value_namer(value)

    def get_file_name_part(self, value):
        return self.file_namer(value)

_PARAMETERS = []

_PARAMETERS.append(_Parameter("quant_method", "Method"))

_PARAMETERS.append(_Parameter(
    "read_depth", "Read depth", is_numeric=True,
    value_namer=lambda x: "{d}x".format(d=x)))

_PARAMETERS.append(_Parameter(
    "read_length", "Read length", is_numeric=True,
    value_namer=lambda x: "{l}b".format(l=x)))

_PARAMETERS.append(_Parameter(
    "paired_end", "End type",
    value_namer=lambda x: "paired-end" if x else "single-end",
    file_namer=lambda x: "pe" if x else "se"))

_PARAMETERS.append(_Parameter(
    "errors", "Error type",
    value_namer=lambda x: "with errors" if x else "no errors",
    file_namer=lambda x: "errors" if x else "no_errors"))

_PARAMETERS.append(_Parameter(
    "bias", "Bias",
    value_namer=lambda x: "with bias" if x else "no bias",
    file_namer=lambda x: "bias" if x else "no_bias"))


def get_parameters():
    return set(_PARAMETERS)


def get_file_name(**params):
    elements = []
    for param in _PARAMETERS:
        if param.name in params:
            value = params[param.name]
            elements.append(param.get_file_name_part(value))
    return "_".join(elements)
