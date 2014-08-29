"""
Utility functions for validating various types of command line options using
the 'schema' package. Exports:

validate_file_option: Check if a file option exists or not.
validate_dir_option: Check if a directory option exists or not.
validate_list_option: Check if a string option is an item in a list
validate_dict_option: Check if a string option is a dictionary key.
validate_options_list: Check if each of a list of items is valid.
validate_int_option: Check if a string option represents an integer.
check_boolean_value: Validates an option string represents a boolean value.
get_logger_for_options: Return a Logger with option-specified severity level.
validate_log_level: Check an option-specified logging level is valid.
substitute_common_options_into_usage: Substitute interpolations into a usage message.
"""

from schema import And, Or, Schema, Use

import log
import os.path
import sys


def validate_file_option(file_option, msg, should_exist=True, nullable=False):
    """
    Check if a file specified by a command line option exists or not.

    Check if the file specified by a command line option exists or not, as
    indicated by the parameter 'should_exist'. The option can be allowed to
    equal 'None' if 'nullable' is set to True. If the specified file fails the
    test, a SchemaError is raised.

    file_option: A string, the path to the file.
    msg: Text for the SchemaError exception raised if the test fails.
    should_exist: Determines if the file is checked for existence or
    non-existence.
    """
    msg = "{msg}: '{file}'.".format(msg=msg, file=file_option)
    validator = open if should_exist else \
        lambda f: not os.path.exists(f)
    if nullable:
        validator = _nullable_validator(validator)
    Schema(validator, error=msg).validate(file_option)


def validate_dir_option(dir_option, msg, should_exist=True, nullable=False):
    """
    Check if a directory specified by a command line option exists or not.

    Check if the directory specified by a command line option exists or not, as
    indicated by the parameter 'should_exist'. The option can be allowed to
    equal 'None' if 'nullable' is set to True. If the specified directory
    fails the test, a SchemaError is raised.

    dir_option: A string, the path to the directory.
    msg: Text for the SchemaError exception raised if the test fails.
    should_exist: Determines if the directory is checked for existence of
    non-existence.
    nullable: If set to True, the command line option is allowed to be 'None'
    (i.e. the option has not been specified).
    """
    msg = "{msg}: '{dir}'.".format(msg=msg, dir=dir_option)
    validator = os.path.isdir if should_exist else \
        lambda f: not os.path.exists(f)
    if nullable:
        validator = _nullable_validator(validator)
    Schema(validator, error=msg).validate(dir_option)


def validate_options_list(option, item_validator, option_name, separator=','):
    """
    Check if each of a list of items is valid according to some validator.

    Check if each item of a command line option (which takes the form of string
    representations of items separated by a specified separator string) is
    valid according to the supplied item validator. Returns a list of objects,
    the items transformed by the validator. If any item is not valid according
    to the validator (i.e. the validator raised an exception), a SchemaError is
    raised.

    option: The command line option, a string.
    item_validator: A type or callable to validate individual items in the
    option string. Should raise an exception if an item is not valid. May
    return transformed versions of items.
    option_name: Name of the option being tested.
    separator: String which separates the command line option string into
    items.
    """
    items = option.split(separator)

    def validated_value(x):
        msg = "'" + str(x) + "' is not a valid " + option_name + "."
        validated = Schema(Use(item_validator), error=msg).validate(x)
        return validated if validated is not None else x

    return [validated_value(i) for i in items]


def validate_list_option(list_option, values_list, msg):
    """
    Check if a command line option is an item in a list.

    Check if a command line option is an item in the specified list. If the
    option is not in the list, a SchemaError is raised.

    list_option: The command line option, a string.
    values_list: A list in which a valid command line option should be an item.
    msg: Text for the SchemaError exception raised if the test fails.
    """
    msg = "{msg}: '{opt}'.".format(msg=msg, opt=list_option)
    Schema(lambda x: x in values_list, error=msg).validate(list_option)


def validate_dict_option(dict_option, values_dict, msg):
    """
    Check if a command line option is a dictionary key.

    Check if a command line option is a key in the specified dictionary and, if
    so, return the corresponding dictionary value. If the option is not a key
    in the dictionary, a SchemaError is raised.

    dict_option: The command line option, a string.
    values_dict: A dictionary, in which a valid command line option should be a
    key.
    msg: Text for the SchemaError exception raised if the test fails.
    """
    msg = "{msg}: '{opt}'.".format(msg=msg, opt=dict_option)
    return Schema(Use(lambda x: values_dict[x]), error=msg).\
        validate(dict_option)


def validate_int_option(int_option, msg, nonneg=False, nullable=False):
    """
    Check if a command line option is an integer.

    Check if a command line option string represents a valid integer and, if
    so, return the integer value. If 'nonneg' is True, the integer must be
    greater than or equal to zero. The option can be allowed to equal 'None' if
    'nullable' is set to True. If the option is not a valid integer, a
    SchemaError is raised.

    int_option: The command line option, a string.
    msg: Text for the SchemaError exception raised if the test fails.
    nonneg: If set to True, the integer must be positive or zero.
    nullable: If set to True, the command line option is allowed to be 'None'
    (i.e. the option has not been specified).
    """
    msg = "{msg}: '{val}'".format(msg=msg, val=int_option)
    validator = Use(int)
    if nonneg:
        validator = And(validator, lambda x: x >= 0)
    if nullable:
        validator = _nullable_validator(validator)

    return Schema(validator, error=msg).validate(int_option)


def check_boolean_value(option_string):
    """
    Validates that a command line option string represents a boolean value.

    Check if a command line option string represents a valid boolean value.
    Accepted strings for True are "true", "t", "yes", "y" or any cased variants
    thereof. Accepted string for False are "false", "f", "no", "n" or any cased
    variants thereof. Any other values are considered invalid, and supplying
    them will cause an exception to be raised.

    option_string: A command line option string representing a boolean value.
    """
    option_string = option_string.lower()
    if option_string in ["true", "t", "yes", "y"]:
        return True
    elif option_string in ["false", "f", "no", "n"]:
        return False
    else:
        raise Exception("Can't convert '{o}' to bool.".format(o=option_string))


def _nullable_validator(validator):
    return Or(validator, None)


def get_logger_for_options(options):
    """
    Return a Logger instance with a command-line specified severity threshold.

    Return a Logger instance with severity threshold specified by the command
    line option log.LOG_LEVEL. Log messages will be written to standard out.
    options: Dictionary mapping from command-line option strings to option
    values.
    """
    return log.get_logger(sys.stderr, options[log.LOG_LEVEL])


def validate_log_level(options):
    """
    Check a command-line option specified logging level is valid.

    Check that the logging level specified through the command-line option
    log.LOG_LEVEL is one of the keys of the LEVELS dictionary.
    options: Dictionary mapping from command-line option strings to option
    values.
    """
    validate_dict_option(
        options[log.LOG_LEVEL], log.LEVELS, "Invalid log level")


def substitute_common_options_into_usage(usage_msg, **substitutions):
    """
    Substitute common option and other interpolations into a usage message.

    Substitute help, version and logging level option specifications and
    descriptions into a script's usage message; also substitute any other
    interpolations specified by additional keyword arguments.
    usage_msg: A script's usage message.
    substitutions: Additional key=value interpolations to substitute into the
    usage message.
    """
    help_spec = "-h --help"
    help_desc = "Show this message."

    ver_spec = "-v --version"
    ver_desc = "Show version."

    log_spec = "--{log_option}=<{log_option}>".format(log_option=log.LOG_LEVEL)
    log_desc = ("Set logging level " +
                "(one of {log_level_vals}) [default: info].").format(
        log_level_vals=str(log.LEVELS.keys()))

    return usage_msg.format(
        log_option_spec=log_spec, log_option_description=log_desc,
        help_option_spec=help_spec, help_option_description=help_desc,
        ver_option_spec=ver_spec, ver_option_description=ver_desc,
        **substitutions)
