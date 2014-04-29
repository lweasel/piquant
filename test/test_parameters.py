import piquant.parameters as parameters


def _get_test_parameter(
        name="name", title="The Name", is_numeric=False, value_namer=None):
    return parameters._Parameter(name, title, is_numeric, value_namer)


def test_get_parameters_returns_parameters_instances():
    params = parameters.get_parameters()
    assert all([isinstance(p, parameters._Parameter) for p in params])


def test_parameter_name_is_correct():
    name = "parameter name"
    p = _get_test_parameter(name=name)
    assert p.name == name


def test_parameter_title_is_correct():
    title = "A Title"
    p = _get_test_parameter(title=title)
    assert p.title == title


def test_parameter_is_numeric_is_correct():
    is_numeric = True
    p = _get_test_parameter(is_numeric=is_numeric)
    assert p.is_numeric == is_numeric


def test_get_value_name_returns_correct_value_when_no_value_namer_supplied():
    value = 30
    p = _get_test_parameter()
    assert p.get_value_name(value) == value


def test_get_value_name_returns_correct_value_when_value_namer_supplied():
    value = 30
    p = _get_test_parameter(value_namer=lambda x: "VAL{v}".format(v=x))
    assert p.get_value_name(value) == "VAL" + str(value)
