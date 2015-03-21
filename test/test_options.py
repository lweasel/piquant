from piquant.options import \
    validate_file_option, validate_dir_option, \
    validate_dict_option, validate_int_option, \
    validate_float_option, validate_options_list, \
    check_boolean_value, validate_list_option
from tempfile import NamedTemporaryFile
from schema import SchemaError
from utils import temp_dir_created

import pytest


def _closed_file_name():
    file = NamedTemporaryFile()
    file_name = file.name
    file.close()
    return file_name


def test_validate_file_option_does_not_raise_exception_for_existing_file_if_file_should_exist():
    with NamedTemporaryFile() as f:
        validate_file_option(f.name, "dummy")


def test_validate_file_option_does_not_raise_exception_for_non_existing_file_if_file_should_not_exist():
    validate_file_option(_closed_file_name(), "dummy", should_exist=False)


def test_validate_file_option_raises_exception_for_non_existing_file_if_file_should_exist():
    with pytest.raises(SchemaError):
        validate_file_option(_closed_file_name(), "dummy")


def test_validate_file_option_raises_exception_for_existing_file_if_file_should_not_exist():
    with NamedTemporaryFile() as f:
        with pytest.raises(SchemaError):
            validate_file_option(f.name, "dummy", should_exist=False)


def test_validate_file_option_raises_exception_for_none_if_nullable_not_specified():
    with pytest.raises(SchemaError):
        validate_file_option(None, "dummy")


def test_validate_file_option_does_not_raise_exception_for_none_if_nullable_specified():
    validate_file_option(None, "dummy", nullable=True)


def test_validate_file_option_exception_message_contains_correct_info():
    file_name = _closed_file_name()

    msg = "dummy"
    with pytest.raises(SchemaError) as exc_info:
        validate_file_option(file_name, msg)

    check_exception_message(exc_info, msg, file_name)


def test_validate_dir_option_does_not_raise_exception_for_existing_dir_if_dir_should_exist():
    with temp_dir_created() as dir_path:
        validate_dir_option(dir_path, "dummy")


def test_validate_dir_option_does_not_raise_exception_for_non_existing_dir_if_dir_should_not_exist():
    dir_path = None
    with temp_dir_created() as dp:
        dir_path = dp

    validate_dir_option(dir_path, "dummy", should_exist=False)


def test_validate_dir_option_raises_exception_for_non_existing_dir_if_dir_should_exist():
    dir_path = None
    with temp_dir_created() as dp:
        dir_path = dp

    with pytest.raises(SchemaError):
        validate_dir_option(dir_path, "dummy")


def test_validate_dir_option_raises_exception_for_existing_dir_if_dir_should_not_exist():
    with temp_dir_created() as dir_path:
        with pytest.raises(SchemaError):
            validate_dir_option(dir_path, "dummy", should_exist=False)


def test_validate_dir_option_raises_exception_for_none_if_nullable_not_specified():
    with pytest.raises(SchemaError):
        validate_dir_option(None, "dummy")


def test_validate_dir_option_does_not_raise_exception_for_none_if_nullable_specified():
    validate_dir_option(None, "dummy", nullable=True)


def test_validate_dir_option_exception_message_contains_correct_info():
    dir_path = None
    with temp_dir_created() as dp:
        dir_path = dp

    msg = "dummy"
    with pytest.raises(SchemaError) as exc_info:
        validate_dir_option(dir_path, msg)

    check_exception_message(exc_info, msg, dir_path)


def test_validate_list_option_does_not_raise_exception_if_item_in_list():
    values_list = [1, 2, 3]
    validate_list_option(values_list[0], values_list, "dummy")


def test_validate_list_option_raises_exception_if_item_not_in_list():
    values_list = [1, 2, 3]
    with pytest.raises(SchemaError):
        validate_list_option(4, values_list, "dummy")


def test_validate_list_option_exception_message_contains_correct_info():
    values_list = [1, 2, 3]

    msg = "dummy"
    item = 4
    with pytest.raises(SchemaError) as exc_info:
        validate_list_option(item, values_list, msg)

    check_exception_message(exc_info, msg, item)


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


def test_validate_int_option_raises_exception_for_non_int():
    with pytest.raises(SchemaError):
        validate_int_option("a", "dummy")


def test_validate_int_option_raises_exception_for_negative_if_min_val_specified():
    with pytest.raises(SchemaError):
        validate_int_option(-1, "dummy", min_val=0)


def test_validate_int_option_does_not_raise_exception_for_negative_if_min_val_not_specified():
    validate_int_option(-1, "dummy")


def test_validate_int_option_raises_exception_for_none_if_nullable_not_specified():
    with pytest.raises(SchemaError):
        validate_int_option(None, "dummy")


def test_validate_int_option_does_not_raise_exception_for_none_if_nullable_specified():
    validate_int_option(None, "dummy", nullable=True)


def test_validate_int_option_exception_message_contains_correct_info():
    msg = "dummy"
    str_val = "abcde"
    with pytest.raises(SchemaError) as exc_info:
        validate_int_option(str_val, msg)

    check_exception_message(exc_info, msg, str_val)


def test_validate_float_option_returns_correct_value():
    float_val = 0.3532
    assert validate_float_option(str(float_val), "dummy") == float_val


def test_validate_float_option_raises_exception_for_non_float():
    with pytest.raises(SchemaError):
        validate_float_option("a", "dummy")


def test_validate_float_option_raises_exception_for_negative_if_min_val_specified():
    with pytest.raises(SchemaError):
        validate_int_option(-1.23, "dummy", min_val=0)


def test_validate_float_option_does_not_raise_exception_for_negative_if_min_val_not_specified():
    validate_int_option(-1.23, "dummy")


def test_validate_int_option_returns_correct_value():
    int_val = 1
    assert validate_int_option(str(int_val), "dummy") == int_val


def test_validate_float_option_exception_message_contains_correct_info():
    msg = "dummy"
    str_val = "abcde"
    with pytest.raises(SchemaError) as exc_info:
        validate_float_option(str_val, msg)

    check_exception_message(exc_info, msg, str_val)


def test_validate_options_list_returns_correct_number_of_items():
    num_items = 5
    separator = ";"
    option_string = (("dummy" + separator) * num_items)[:-1]
    assert len(validate_options_list(
        option_string, lambda x: x, "dummy_name", separator)) == num_items


def test_validate_options_list_returns_transformed_objects():
    option_values = [1, 5, 10]
    option_string = ",".join([str(v) for v in option_values])
    assert validate_options_list(
        option_string, int, "dummy_name") == option_values


def test_validate_options_list_raises_exception_for_invalid_value():
    option_values = [1, 5, "ten"]
    option_string = ",".join([str(v) for v in option_values])
    with pytest.raises(SchemaError):
        validate_options_list(option_string, int, "dummy_name")


def test_validate_options_list_exception_message_contains_correct_info():
    option_name = "option_name"
    invalid_option = "Three"
    option_values = [1, 2, invalid_option]
    option_string = ",".join([str(v) for v in option_values])

    with pytest.raises(SchemaError) as exc_info:
        validate_options_list(option_string, int, option_name)

    check_exception_message(exc_info, option_name, invalid_option)


def test_check_boolean_value_accepts_valid_true_strings():
    for option_string in ["true", "t", "yes", "y"]:
        assert check_boolean_value(option_string)
        assert check_boolean_value(option_string.upper())


def test_check_boolean_value_accepts_valid_false_strings():
    for option_string in ["false", "f", "no", "n"]:
        assert not check_boolean_value(option_string)
        assert not check_boolean_value(option_string.upper())


def test_check_boolean_value_raises_exception_for_invalid_string():
    with pytest.raises(Exception):
        check_boolean_value("not a boolean")


def check_exception_message(exc_info, *args):
    exc_msg = str(exc_info)
    for arg in args:
        assert str(arg) in exc_msg
