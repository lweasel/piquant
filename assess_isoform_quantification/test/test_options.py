from assess_isoform_quantification.options \
    import validate_file_option, validate_dict_option, validate_int_option
from tempfile import NamedTemporaryFile
from schema import SchemaError

import pytest


def test_validate_file_option_returns_handle_for_file_that_exists():
    file = NamedTemporaryFile()
    file_name = file.name
    assert validate_file_option(file_name, "dummy") is not None


def test_validate_file_option_raises_exception_for_non_existing_file():
    file = NamedTemporaryFile()
    file_name = file.name
    file.close()
    with pytest.raises(SchemaError):
        validate_file_option(file_name, "dummy")


def test_validate_file_option_exception_message_contains_correct_info():
    file = NamedTemporaryFile()
    file_name = file.name
    file.close()

    msg = "dummy"
    with pytest.raises(SchemaError) as exc_info:
        validate_file_option(file_name, msg)

    check_exception_message(exc_info, msg, file_name)


def test_validate_dict_option_returns_correct_dict_value():
    values_dict = {1: "one", 2: "two"}
    assert validate_dict_option(2, values_dict, "dummy") == "two"


def test_validate_dict_option_raises_exception_for_non_existent_key():
    values_dict = {1: "one", 2: "two"}
    with pytest.raises(SchemaError):
        validate_dict_option(3, values_dict, "dummy")


def test_validate_dict_option_exception_message_contains_correct_info():
    values_dict = {1: "one", 2: "two"}

    msg = "dummy"
    key = 3
    with pytest.raises(SchemaError) as exc_info:
        validate_dict_option(key, values_dict, msg)

    check_exception_message(exc_info, msg, key)


def test_validate_int_option_returns_correct_value():
    int_val = 1
    assert validate_int_option(str(1), "dummy") == int_val


def test_validate_int_option_raises_exception_for_non_int():
    with pytest.raises(SchemaError):
        validate_int_option("a", "dummy")


def test_validate_int_option_raises_exception_for_negative_if_nonneg_specified():
    with pytest.raises(SchemaError):
        validate_int_option(-1, "dummy", nonneg=True)


def test_validate_int_option_does_not_raise_exception_for_negative_if_nonneg_not_specified():
    validate_int_option(-1, "dummy")


def test_validate_int_option_exception_message_contains_correct_info():
    msg = "dummy"
    str_val = "abcde"
    with pytest.raises(SchemaError) as exc_info:
        validate_int_option(str_val, msg)

    check_exception_message(exc_info, msg, str_val)


def check_exception_message(exc_info, *args):
    exc_msg = exc_info.value.message
    for arg in args:
        assert str(arg) in exc_msg
