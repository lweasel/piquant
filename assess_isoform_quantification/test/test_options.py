from assess_isoform_quantification.options \
    import validate_file_option, validate_dict_option
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


def check_exception_message(exc_info, *args):
    exc_msg = exc_info.value.message
    for arg in args:
        assert str(arg) in exc_msg
