from assess_isoform_quantification.options import validate_file_option
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

    exc_msg = exc_info.value.message
    assert msg in exc_msg
    assert file_name in exc_msg
